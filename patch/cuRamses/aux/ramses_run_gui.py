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


class StageEnd(Exception):
    """The next configuration stage has been reached."""


def answer_text(question, value):
    if question.kind == 'edit':
        return json.dumps(value, indent=2)
    if question.kind == 'floats' and isinstance(value, list):
        return ','.join(map(str, value))
    if question.kind == 'bool':
        return value
    return '' if value is None else str(value)


class StageUI(ReplayUI):
    """Collect one whole form through the same generator as the terminal UI.

    Drafts are keyed by question identity, not position: changing a branch
    cannot feed an old answer into a different parameter.
    """
    def __init__(self, generator, answers=(), draft=None):
        super().__init__(answers)
        self.generator = generator
        self.draft = draft or {}
        self.section = 'Run files'
        self.stage = None
        self.fields = []
        self.values = []
        self.error = None

    def info(self, text=''):
        super().info(text)
        heading = str(text).strip()
        if heading.startswith('==='):
            self.section = heading.strip('= ')
            if self.section.startswith('lagRamses run generator'):
                self.section = 'Run files'

    def _answer(self, question):
        if question.prompt.startswith('==='):
            self.section = question.prompt.strip('= ')
        elif question.prompt == 'Zoom-in run?':
            self.section = 'Zoom region'
        elif 'full parameter editor' in question.prompt:
            self.section = 'Advanced settings'
        if self.index < len(self.answers):
            return super()._answer(question)
        if self.stage is not None and self.section != self.stage:
            raise StageEnd()
        self.stage = self.section
        key = (question.kind, question.prompt)
        raw = self.draft.get(key, answer_text(question, question.default))
        self.fields.append((key, question, raw))
        value = parse_answer(question, raw, self.generator)
        self.values.append(value)
        return copy.deepcopy(value)


def collect_stage(generator, answers=(), draft=None):
    ui = StageUI(generator, answers, draft)
    files = OrderedDict()
    report = None
    try:
        report = generator.generate_run(ui, files.__setitem__)
    except StageEnd:
        pass
    except (ValueError, KeyError, TypeError, OverflowError) as exc:
        ui.error = str(exc)
    return ui, files, report


class RunWizard:
    def __init__(self, root, generator, tk, ttk, dialogs):
        self.root, self.generator, self.tk, self.ttk = root, generator, tk, ttk
        self.filedialog, self.messagebox = dialogs
        self.answers, self.files, self.report = [], OrderedDict(), None
        self.history = []
        self.draft = {}
        root.title('lagRamses run setup')
        root.geometry('1100x800')
        root.minsize(760, 560)
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
        if self.history:
            count, self.draft = self.history.pop()
            del self.answers[count:]
        self.refresh()

    def refresh(self):
        for child in self.frame.winfo_children():
            child.destroy()
        self.files, self.report = OrderedDict(), None
        ui, files, report = collect_stage(self.generator, self.answers, self.draft)
        self.stage_ui = ui
        if not ui.fields and not ui.error:
            self.files, self.report = files, report
            self.show_preview(ui.notes)
            return
        self.show_stage(ui)

    def read_draft(self):
        return {key: read() for key, read in self.readers.items()}

    def update_fields(self):
        self.draft = self.read_draft()
        self.refresh()

    def browse(self, variable, directory=False):
        chooser = self.filedialog.askdirectory if directory else self.filedialog.askopenfilename
        path = chooser(parent=self.root, mustexist=True) if directory else chooser(parent=self.root)
        if path:
            variable.set(path)

    def show_stage(self, ui):
        self.ttk.Label(self.frame, text='{} — {}'.format(len(self.history) + 1, ui.stage or 'Setup'),
                       font=('', 16)).pack(anchor='w', pady=(0, 12))
        self.ttk.Label(self.frame, text='Edit the settings together. Model choices update related fields.',
                       wraplength=950).pack(anchor='w')

        # Keep navigation visible while long model menus/advanced edits scroll.
        bar = self.ttk.Frame(self.frame)
        bar.pack(side='bottom', fill='x', pady=10)
        self.ttk.Button(bar, text='Back', command=self.back,
                        state='normal' if self.history else 'disabled').pack(side='left')
        self.ttk.Button(bar, text='Next', command=self.next).pack(side='right')
        self.ttk.Button(bar, text='Update fields', command=self.update_fields).pack(side='right', padx=8)
        self.error = self.ttk.Label(self.frame, text=ui.error or '', foreground='firebrick', wraplength=950)
        self.error.pack(side='bottom', anchor='w', pady=8)

        holder = self.ttk.Frame(self.frame)
        holder.pack(fill='both', expand=True)
        canvas = self.tk.Canvas(holder, highlightthickness=0)
        scrollbar = self.ttk.Scrollbar(holder, orient='vertical', command=canvas.yview)
        canvas.configure(yscrollcommand=scrollbar.set)
        scrollbar.pack(side='right', fill='y')
        canvas.pack(side='left', fill='both', expand=True)
        content = self.ttk.Frame(canvas)
        window = canvas.create_window((0, 0), window=content, anchor='nw')
        content.bind('<Configure>', lambda event: canvas.configure(scrollregion=canvas.bbox('all')))
        canvas.bind('<Configure>', lambda event: canvas.itemconfigure(window, width=event.width))
        for column in range(2):
            content.columnconfigure(column, weight=1, uniform='fields')
        self.readers, self.controls = {}, {}
        cell = 0
        for key, question, raw in ui.fields:
            wide = question.kind in ('choice', 'edit')
            if wide and cell % 2:
                cell += 1
            panel = self.ttk.LabelFrame(content, text=question.prompt.strip('= \n'), padding=10)
            panel.grid(row=cell // 2, column=0 if wide else cell % 2,
                       columnspan=2 if wide else 1, sticky='nsew', padx=5, pady=5)
            cell += 2 if wide else 1
            if question.kind == 'edit':
                entry = self.tk.Text(panel, height=16, wrap='none')
                entry.pack(fill='both', expand=True)
                entry.insert('1.0', raw)
                self.readers[key] = lambda entry=entry: entry.get('1.0', 'end-1c')
                self.controls[key] = entry
                continue
            variable = (self.tk.BooleanVar(value=raw) if question.kind == 'bool'
                        else self.tk.StringVar(value=raw))
            self.controls[key] = variable
            self.readers[key] = variable.get
            if question.kind == 'choice':
                for column in range(2):
                    panel.columnconfigure(column, weight=1, uniform='choices')
                for index, (value, labels) in enumerate(question.options.items()):
                    self.ttk.Radiobutton(panel, text=labels[0], value=value, variable=variable,
                                         command=self.update_fields).grid(
                                             row=index // 2, column=index % 2, sticky='w', padx=5, pady=5)
            elif question.kind == 'bool':
                for index, (label, value) in enumerate((('Yes', True), ('No', False))):
                    self.ttk.Radiobutton(panel, text=label, value=value, variable=variable,
                                         command=self.update_fields).grid(row=0, column=index, padx=12)
            else:
                entry = self.ttk.Entry(panel, textvariable=variable)
                entry.pack(fill='x')
                if question.prompt == 'Output directory' or 'path' in question.prompt.lower():
                    self.ttk.Button(panel, text='Browse…', command=lambda var=variable,
                                    directory=question.prompt == 'Output directory':
                                    self.browse(var, directory)).pack(anchor='w', pady=5)

    def next(self):
        draft = self.read_draft()
        ui, _, _ = collect_stage(self.generator, self.answers, draft)
        if ui.error:
            self.error.configure(text=ui.error)
            return
        # A value may reveal a dependent field. Show it for review before
        # accepting newly exposed defaults (e.g. quintessence potential).
        if [key for key, _, _ in ui.fields] != list(self.readers):
            self.draft = draft
            self.refresh()
            return
        self.history.append((len(self.answers), draft))
        self.answers.extend(ui.values)
        self.draft = {}
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
              'Tk support, or run mkrun.py --mode cli for the CLI.'.format(exc), file=sys.stderr)
        return 2
    try:
        root = tk.Tk()
    except tk.TclError as exc:
        print('GUI unavailable: cannot open a graphical display ({}). Run from '
              'a desktop/display-enabled session, or use mkrun.py --mode cli.'.format(exc),
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
