#!/usr/bin/env python3
"""Tkinter frontend for mkrun's shared wizard; imported only for --gui.

The wizard is replayed from recorded answers after each step. Its output sink
stores strings in memory, so Back, browsing and preview never write files.
No tkinter import or display is required for the model and save tests.
"""
import copy
import json
import math
import os
import re
import stat
import tempfile
from collections import OrderedDict


class Question(Exception):
    def __init__(self, kind, prompt, default, cast=str, options=None):
        super().__init__(prompt)
        self.kind, self.prompt, self.default = kind, prompt, default
        self.cast, self.options = cast, options


class ReplayUI:
    """Injectable prompts for the exact same CLI collection/generation path."""
    def __init__(self, answers=()):
        self.answers = copy.deepcopy(list(answers))
        self.index = 0
        self.notes = []

    def _answer(self, question):
        if self.index == len(self.answers):
            raise question
        value = self.answers[self.index]
        self.index += 1
        return copy.deepcopy(value)

    def ask(self, prompt, default=None, cast=str):
        return self._answer(Question('value', prompt, default, cast))

    def ask_bool(self, prompt, default=True):
        return self._answer(Question('bool', prompt.strip(), default))

    def ask_choice(self, prompt, options, default_key):
        return self._answer(Question('choice', prompt.strip() or 'IC pipeline',
                                     default_key, options=options))

    def ask_floats(self, prompt, default_csv):
        return self._answer(Question('floats', prompt, default_csv))

    def edit(self, values):
        return self._answer(Question('edit', 'Advanced namelist parameters', values))

    def info(self, text=''):
        self.notes.append(str(text))


def parse_answer(question, raw, generator):
    """Strict syntax checks before replaying answers into the shared workflow."""
    if question.kind == 'bool':
        return bool(raw)
    if question.kind == 'choice':
        if raw not in question.options:
            raise ValueError('Select one of the listed options.')
        return raw
    if question.kind == 'edit':
        values = json.loads(raw, object_pairs_hook=OrderedDict)
        if not isinstance(values, dict):
            raise ValueError('Advanced parameters must be a JSON object.')
        result = OrderedDict()
        for key, value in values.items():
            name = key.lower()
            param = generator.rng.PARAM_BY_NAME.get(name)
            if param is None:
                raise ValueError('Unknown parameter: {}'.format(key))
            if name in result:
                raise ValueError('Duplicate parameter: {}'.format(key))
            validate_value(param.ftype, value, key)
            result[name] = value
        return result
    raw = raw.strip()
    if not raw:
        if question.default is None:
            raise ValueError('A value is required.')
        raw = str(question.default)
    if question.kind == 'floats':
        values = [float(token.strip()) for token in raw.split(',')]
        if not all(math.isfinite(value) for value in values):
            raise ValueError('Redshifts must be finite numbers.')
        return values
    value = question.cast(raw)
    if isinstance(value, float) and not math.isfinite(value):
        raise ValueError('Enter a finite number.')
    if isinstance(value, str) and any(c in value for c in '\n\r\x00'):
        raise ValueError('Enter a single line without NUL characters.')
    # The shared legacy writers interpolate strings into quoted Fortran/INI.
    if isinstance(value, str) and any(c in value for c in "'\"!#"):
        raise ValueError('Quotes and comment markers are not supported in this field.')
    return value


def validate_value(ftype, value, name):
    valid = True
    if ftype == 'bool':
        valid = type(value) is bool
    elif ftype == 'int':
        valid = type(value) is int
    elif ftype == 'real':
        valid = type(value) in (int, float) and math.isfinite(value)
    elif ftype in ('int_arr', 'real_arr'):
        valid = isinstance(value, str)
        if valid and value.strip():
            scalar = r'[+-]?\d+' if ftype == 'int_arr' else (
                r'[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eEdD][+-]?\d+)?')
            valid = all(re.fullmatch(r'(?:[1-9]\d*\*)?' + scalar, part.strip())
                        for part in value.split(','))
            if valid and ftype == 'real_arr':
                valid = all(math.isfinite(float(part.split('*')[-1].lower().replace('d', 'e')))
                            for part in value.split(','))
    else:
        valid = isinstance(value, str) and not any(c in value for c in '\n\r\x00')
        if valid and ftype == 'str':
            # Existing CLI solver choices may already be single quoted.
            body = value[1:-1] if value.startswith("'") and value.endswith("'") else value
            valid = not any(c in body for c in "'\"!#")
    if not valid:
        raise ValueError('{} must be a valid {} value.'.format(name, ftype))


def snapshot_targets(files):
    """Read only: capture exact config bytes for confirmation and race checks."""
    result = {}
    for path in files:
        if not os.path.lexists(path):
            result[path] = None
            continue
        mode = os.lstat(path)
        if not stat.S_ISREG(mode.st_mode) or mode.st_size > 4 * 1024 * 1024:
            raise ValueError('Refusing to replace a symlink, non-file or large file: ' + path)
        with open(path, 'rb') as stream:
            result[path] = stream.read()
    return result


def save_preview(files, expected):
    """Save exactly the confirmed preview, refusing targets changed meanwhile.

    Replacements are atomic per file. On an I/O failure the error names any
    files already saved; a multi-file bundle is not a filesystem transaction.
    """
    if not files or set(files) != set(expected):
        raise ValueError('No complete confirmed preview to save.')
    if snapshot_targets(files) != expected:
        raise ValueError('A destination changed. Review and confirm the preview again.')
    saved = []
    try:
        for path, text in files.items():
            directory = os.path.dirname(path)
            os.makedirs(directory, exist_ok=True)
            fd, temp_path = tempfile.mkstemp(prefix='.mkrun-', dir=directory)
            try:
                with os.fdopen(fd, 'w', encoding='utf-8') as stream:
                    stream.write(text)
                if expected[path] is None:
                    # Publish only a complete file; link refuses a concurrent creation.
                    os.link(temp_path, path)
                else:
                    if snapshot_targets({path: text})[path] != expected[path]:
                        raise ValueError('Destination changed: ' + path)
                    os.chmod(temp_path, stat.S_IMODE(os.stat(path).st_mode))
                    os.replace(temp_path, path)
            finally:
                if os.path.exists(temp_path):
                    os.unlink(temp_path)
            saved.append(path)
    except (OSError, ValueError) as exc:
        raise ValueError('{}\nFiles already saved: {}'.format(
            exc, ', '.join(saved) or 'none')) from exc
    return saved


class RunWizard:
    def __init__(self, root, generator, tk, ttk, dialogs):
        self.root, self.generator, self.tk, self.ttk = root, generator, tk, ttk
        self.filedialog, self.messagebox = dialogs
        self.answers, self.files, self.report = [], OrderedDict(), None
        self.previous_answer = None
        self.returning = False
        root.title('lagRamses run setup')
        root.geometry('980x720')
        root.minsize(720, 520)
        self.frame = ttk.Frame(root, padding=16)
        self.frame.pack(fill='both', expand=True)
        root.protocol('WM_DELETE_WINDOW', self.close)
        self.refresh()

    def close(self):
        if not self.answers or self.messagebox.askyesno(
                'Close setup', 'Close this setup wizard? Unsaved changes will be lost.',
                parent=self.root):
            self.root.destroy()

    def back(self):
        if self.answers:
            self.previous_answer = self.answers.pop()
            self.returning = True
        self.refresh()

    def refresh(self):
        for child in self.frame.winfo_children():
            child.destroy()
        self.files, self.report = OrderedDict(), None
        ui = ReplayUI(self.answers)
        try:
            report = self.generator.generate_run(ui, self.files.__setitem__)
        except Question as question:
            self.show_question(question, ui.notes)
            return
        except (ValueError, KeyError, TypeError, OverflowError) as exc:
            self.ttk.Label(self.frame, text='Setup needs correction', font=('', 16)).pack(anchor='w')
            self.ttk.Label(self.frame, text=str(exc), wraplength=850).pack(anchor='w', pady=20)
            self.ttk.Button(self.frame, text='Back', command=self.back).pack(anchor='w')
            return
        self.report = report
        self.show_preview(ui.notes)

    def show_question(self, question, notes):
        if self.returning:
            question.default = self.previous_answer
            if question.kind == 'floats':
                question.default = ','.join(str(value) for value in question.default)
            self.returning = False
        self.question = question
        self.ttk.Label(self.frame, text='Run setup — step {}'.format(len(self.answers) + 1),
                       font=('', 16)).pack(anchor='w')
        headings = [note.strip() for note in notes if note.strip().startswith('===')]
        if headings:
            self.ttk.Label(self.frame, text=headings[-1].strip('= ')).pack(anchor='w', pady=8)
        self.ttk.Label(self.frame, text=question.prompt.strip('= \n'),
                       wraplength=850).pack(anchor='w', pady=12)
        if question.kind == 'edit':
            self.ttk.Label(self.frame, text='Edit named parameters as JSON. Booleans use true/false; '
                           'arrays use quoted Fortran lists (for example "8*8.").\n'
                           'The shared wizard subsequently sets its output schedule, IC paths '
                           'and selected cooling defaults.', wraplength=850).pack(anchor='w')
            self.entry = self.tk.Text(self.frame, height=16, wrap='none')
            self.entry.pack(fill='both', expand=True, pady=8)
            self.entry.insert('1.0', json.dumps(question.default, indent=2))
        elif question.kind == 'bool':
            self.variable = self.tk.BooleanVar(value=question.default)
            self.ttk.Checkbutton(self.frame, text='Yes', variable=self.variable).pack(anchor='w')
        elif question.kind == 'choice':
            self.variable = self.tk.StringVar(value=question.default)
            # A scrollable list keeps the gravity menu usable on small displays.
            holder = self.ttk.Frame(self.frame)
            holder.pack(fill='both', expand=True)
            self.choice_list = self.tk.Listbox(holder, exportselection=False, height=10)
            scrollbar = self.ttk.Scrollbar(holder, command=self.choice_list.yview)
            self.choice_list.configure(yscrollcommand=scrollbar.set)
            scrollbar.pack(side='right', fill='y')
            self.choice_list.pack(side='left', fill='both', expand=True)
            for index, (key, labels) in enumerate(question.options.items()):
                self.choice_list.insert('end', '{} — {}'.format(key, labels[0]))
                if key == question.default:
                    self.choice_list.selection_set(index)
                    self.choice_list.see(index)
        else:
            self.variable = self.tk.StringVar(value='' if question.default is None else str(question.default))
            self.entry = self.ttk.Entry(self.frame, textvariable=self.variable)
            self.entry.pack(fill='x', pady=8)
            self.entry.focus_set()
            self.entry.bind('<Return>', lambda event: self.next())
            if question.prompt == 'Output directory':
                self.ttk.Button(self.frame, text='Choose existing directory…',
                                command=self.browse_directory).pack(anchor='w')
            elif 'path' in question.prompt.lower():
                self.ttk.Button(self.frame, text='Choose existing file…',
                                command=self.browse_file).pack(anchor='w')
        self.error = self.ttk.Label(self.frame, text='', foreground='firebrick', wraplength=850)
        self.error.pack(anchor='w', pady=12)
        buttons = self.ttk.Frame(self.frame)
        buttons.pack(side='bottom', fill='x', pady=10)
        self.ttk.Button(buttons, text='Back', command=self.back,
                        state='normal' if self.answers else 'disabled').pack(side='left')
        self.ttk.Button(buttons, text='Next', command=self.next).pack(side='right')
        self.ttk.Label(self.frame, text='Configuration files only. Saving is a separate step.').pack(
            side='bottom', anchor='w')

    def browse_directory(self):
        path = self.filedialog.askdirectory(parent=self.root, mustexist=True)
        if path:
            self.variable.set(path)

    def browse_file(self):
        path = self.filedialog.askopenfilename(parent=self.root)
        if path:
            self.variable.set(path)

    def next(self):
        try:
            if self.question.kind == 'edit':
                raw = self.entry.get('1.0', 'end-1c')
            elif self.question.kind == 'choice':
                selection = self.choice_list.curselection()
                if not selection:
                    raise ValueError('Select an option.')
                raw = list(self.question.options)[selection[0]]
            else:
                raw = self.variable.get()
            answer = parse_answer(self.question, raw, self.generator)
        except (ValueError, TypeError) as exc:
            self.error.configure(text=str(exc))
            return
        self.answers.append(answer)
        self.refresh()

    def show_preview(self, notes):
        self.ttk.Label(self.frame, text='Review generated configuration', font=('', 16)).pack(anchor='w')
        self.ttk.Label(self.frame, text=self.report['outdir'], wraplength=850).pack(anchor='w', pady=8)
        book = self.ttk.Notebook(self.frame)
        book.pack(fill='both', expand=True)
        validation = '\n'.join(str(msg) for msg in self.report['messages']) or 'No errors from the shared parameter validator.'
        summary = (validation + '\n\nConfiguration preview only; this is not a launch approval.\n'
                   'Review the output schedule and storage requirements before running.\n\n' +
                   '\n'.join(notes))
        for label, text in [('Validation', summary)] + [
                (os.path.basename(path), text) for path, text in self.files.items()]:
            panel = self.ttk.Frame(book)
            book.add(panel, text=label)
            widget = self.tk.Text(panel, wrap='none')
            vertical = self.ttk.Scrollbar(panel, orient='vertical', command=widget.yview)
            horizontal = self.ttk.Scrollbar(panel, orient='horizontal', command=widget.xview)
            widget.configure(yscrollcommand=vertical.set, xscrollcommand=horizontal.set)
            vertical.pack(side='right', fill='y')
            horizontal.pack(side='bottom', fill='x')
            widget.pack(fill='both', expand=True)
            widget.insert('1.0', text)
            widget.configure(state='disabled')
        bar = self.ttk.Frame(self.frame)
        bar.pack(fill='x', pady=12)
        self.ttk.Button(bar, text='Back', command=self.back).pack(side='left')
        errors = any(msg.level == 'ERROR' for msg in self.report['messages'])
        self.save_button = self.ttk.Button(bar, text='Save files…', command=self.save,
                                           state='disabled' if errors else 'normal')
        self.save_button.pack(side='right')

    def save(self):
        if not self.report or any(msg.level == 'ERROR' for msg in self.report['messages']):
            return
        try:
            expected = snapshot_targets(self.files)
            existing = [path for path, content in expected.items() if content is not None]
            prompt = 'Save these configuration files?\n\n' + '\n'.join(self.files)
            if existing:
                prompt += '\n\nOVERWRITE existing files:\n' + '\n'.join(existing)
            if not self.messagebox.askyesno('Confirm file writes', prompt, parent=self.root):
                return
            saved = save_preview(self.files, expected)
        except (OSError, ValueError) as exc:
            self.messagebox.showerror('Cannot save configuration', str(exc), parent=self.root)
            return
        self.save_button.configure(state='disabled')
        self.messagebox.showinfo('Configuration saved', '\n'.join(saved), parent=self.root)


def launch(generator):
    """Return a normal CLI status even when Tk or a graphical display is absent."""
    import sys
    try:
        import tkinter as tk
        from tkinter import filedialog, messagebox, ttk
    except ImportError as exc:
        print('GUI unavailable: Tkinter is not installed ({}). Use a Python with '
              'Tk support, or run mkrun.py without --gui for the CLI.'.format(exc), file=sys.stderr)
        return 2
    try:
        root = tk.Tk()
    except tk.TclError as exc:
        print('GUI unavailable: cannot open a graphical display ({}). Run from '
              'a desktop/display-enabled session, or omit --gui for the CLI.'.format(exc),
              file=sys.stderr)
        return 2
    try:
        RunWizard(root, generator, tk, ttk, (filedialog, messagebox))
        root.mainloop()
    finally:
        try:
            root.destroy()
        except tk.TclError:
            pass
    return 0
