#!/usr/bin/env python3
"""Bounded setup-only tests: python3 -B -m unittest discover -s patch/cuRamses/aux -p test_ramses_run_gui.py -v."""
import contextlib
import io
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import types
import unittest
from unittest import mock

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
import mkrun
import ramses_run_gui as gui


def collect(overrides=None):
    """Walk every real shared prompt, using its own defaults unless overridden."""
    overrides = overrides or {}
    answers = []
    while True:
        ui = gui.ReplayUI(answers)
        files = {}
        try:
            report = mkrun.generate_run(ui, files.__setitem__)
            return answers, files, report
        except gui.Question as question:
            key = question.prompt.strip('\n= ')
            value = overrides.get(key, question.default)
            if question.kind == 'floats' and isinstance(value, str):
                value = gui.parse_answer(question, value, mkrun)
            answers.append(value)


class WizardTests(unittest.TestCase):
    def test_sink_formation_prerequisites_and_round_trip(self):
        rng = mkrun.rng
        for hydro in (False, True):
            for poisson in (False, True):
                for create in (None, False, True):
                    values = dict(sink=True, hydro=hydro, poisson=poisson)
                    if create is not None:
                        values['create_sinks'] = create
                    with self.subTest(**values):
                        errors = [msg for msg in rng.validate_params(values)
                                  if msg.level == 'ERROR']
                        self.assertEqual(bool(errors),
                                         create is not False and not (hydro and poisson))
        values = dict(sink=True, hydro=True, poisson=False, create_sinks=False)
        text = rng.format_namelist(values)
        self.assertIn('&SINK_PARAMS', text)
        parsed, _ = rng.parse_namelist(text)
        self.assertIs(rng.import_to_values(parsed)['create_sinks'], False)
        self.assertNotIn('&SINK_PARAMS', rng.format_namelist(dict(sink=False)))

    def test_invalid_sink_edit_does_not_write(self):
        answers, _, _ = collect({
            'Open the full parameter editor for fine-tuning before writing?': True})
        ui = gui.ReplayUI(answers)
        with mock.patch.object(ui, 'edit', return_value={
                'sink': True, 'hydro': True, 'poisson': False}), \
                mock.patch('builtins.open', side_effect=AssertionError('unexpected write')):
            files = {}
            with self.assertRaisesRegex(ValueError, 'sink creation requires'):
                mkrun.generate_run(ui, files.__setitem__)
            self.assertEqual(files, {})

    def test_preview_never_touches_filesystem(self):
        with mock.patch('builtins.open', side_effect=AssertionError('preview wrote a file')), \
                mock.patch('os.makedirs', side_effect=AssertionError('preview made a directory')):
            answers, files, report = collect()
        self.assertEqual(len(files), 3)
        self.assertTrue(report['values']['cosmo'])
        self.assertFalse(report['values']['hydro'])
        self.assertTrue(answers)

    def test_cli_and_gui_same_bytes_for_all_ic_pipelines(self):
        cases = [(False, 'music'), (False, 'monofonic'), (False, 'none'),
                 (True, 'music'), (True, 'genetic'), (True, 'genetic_mono'), (True, 'none')]
        for hydro in (False, True):
            for zoom, pipeline in cases:
                with self.subTest(hydro=hydro, zoom=zoom, pipeline=pipeline):
                    overrides = {'Run mode': 'hydro' if hydro else 'dmo',
                                 'Zoom-in run?': zoom, 'IC pipeline': pipeline}
                    answers, preview, _ = collect(overrides)
                    response = iter(answers)

                    def terminal_input(prompt):
                        value = next(response)
                        if isinstance(value, bool):
                            return 'yes' if value else 'no'
                        if isinstance(value, list):
                            return ','.join(map(str, value))
                        return str(value)

                    cli = {}
                    with mock.patch('builtins.input', side_effect=terminal_input), \
                            contextlib.redirect_stdout(io.StringIO()):
                        mkrun.generate_run(write_text=cli.__setitem__)
                    self.assertEqual(preview, cli)

    def test_valid_outputs_match_head_before_refactor(self):
        # The baseline is loaded into memory, never restored over the worktree.
        baseline = 'b1d489633822c4ecca2cd9c68cc5b592b4ec25f6:mkrun.py'
        source = subprocess.check_output(['git', 'show', baseline], cwd=str(ROOT), text=True)
        old = types.ModuleType('mkrun_head_reference')
        old.__file__ = str(ROOT / 'mkrun.py')
        exec(compile(source, 'HEAD:mkrun.py', 'exec'), old.__dict__)
        for overrides in ({}, {'Run mode': 'hydro'},
                          {'Zoom-in run?': True, 'IC pipeline': 'genetic_mono'},
                          {'Dark matter sector': 'sidm', 'Gravity / dark-energy sector': 'fR'}):
            with self.subTest(overrides=overrides), tempfile.TemporaryDirectory() as directory:
                settings = dict(overrides, **{'Output directory': directory})
                answers, preview, _ = collect(settings)
                response = iter(answers)

                def terminal_input(prompt):
                    value = next(response)
                    if isinstance(value, bool):
                        return 'yes' if value else 'no'
                    if isinstance(value, list):
                        return ','.join(map(str, value))
                    return str(value)

                with mock.patch('builtins.input', side_effect=terminal_input), \
                        contextlib.redirect_stdout(io.StringIO()):
                    old.main()
                self.assertEqual(preview, {path: Path(path).read_text() for path in preview})

    def test_invalid_inputs_do_not_write(self):
        for settings in ({'Run name (used as file/dir prefix)': '../escape'},
                         {'box size [Mpc/h]': 0},
                         {'output redshifts (comma separated, high-z first)': '-1'},
                         {'output redshifts (comma separated, high-z first)': '100'},
                         {'levelmin (base/coarse level)': 31},
                         {'omega_b': 0.9}):
            with self.subTest(settings=settings), \
                    mock.patch('builtins.open', side_effect=AssertionError('unexpected write')):
                with self.assertRaises(ValueError):
                    collect(settings)

    def test_missing_database_sector_is_clear_error(self):
        if 'm_axion' in mkrun.rng.PARAM_BY_NAME:
            self.skipTest('Current database supports FDM')
        with self.assertRaisesRegex(ValueError, 'does not support m_axion'):
            collect({'Dark matter sector': 'fdm'})

    def test_strict_numeric_and_advanced_inputs(self):
        q = gui.Question('value', 'number', 1.0, float)
        for raw in ('nan', 'inf', 'not-a-number'):
            with self.assertRaises(ValueError):
                gui.parse_answer(q, raw, mkrun)
        q = gui.Question('edit', 'advanced', {})
        for raw in ('[]', '{"typo":1}', '{"hydro":"false"}',
                    '{"levelmin":1.5}', '{"aout":"nan"}', '{"aout":"1,2/ &RUN_PARAMS"}'):
            with self.assertRaises(ValueError):
                gui.parse_answer(q, raw, mkrun)
        values = gui.parse_answer(q, '{"hydro":false,"levelmin":8,"m_refine":"8*8."}', mkrun)
        self.assertIs(values['hydro'], False)

    def test_back_replay_accepts_new_branch_without_stale_files(self):
        answers, _, _ = collect()
        # First two answers are name/output directory; third is run mode.
        with self.assertRaises(gui.Question) as caught:
            mkrun.generate_run(gui.ReplayUI(answers[:2]), {}.__setitem__)
        self.assertEqual(caught.exception.default, 'dmo')
        _, files, report = collect({'Run mode': 'hydro', 'IC pipeline': 'none'})
        self.assertEqual(len(files), 1)
        self.assertTrue(report['values']['hydro'])

    def test_help_needs_no_tk(self):
        with mock.patch.dict(sys.modules, {'tkinter': None}), contextlib.redirect_stdout(io.StringIO()):
            with self.assertRaises(SystemExit) as caught:
                mkrun.main(['--help'])
        self.assertEqual(caught.exception.code, 0)

    def test_mode_gui_dispatches_to_shared_frontend(self):
        with mock.patch('ramses_run_gui.launch', return_value=0) as launch:
            self.assertEqual(mkrun.main(['--mode', 'gui']), 0)
            self.assertEqual(mkrun.main(['--gui', '--mode', 'gui']), 0)
        self.assertEqual(launch.call_count, 2)
        self.assertTrue(all(call == mock.call(mkrun) for call in launch.call_args_list))

    def test_grouped_forms_match_shared_generator(self):
        for overrides in ({}, {'Run mode': 'hydro', 'Zoom-in run?': True},
                          {'Dark matter sector': 'sidm', 'Gravity / dark-energy sector': 'w0wa'},
                          {'Zoom-in run?': True, 'IC pipeline': 'genetic_mono'}):
            with self.subTest(overrides=overrides):
                expected_answers, expected_files, _ = collect(overrides)
                answers, stages = [], {}
                for _ in range(20):
                    ui, _, _ = gui.collect_stage(mkrun, answers)
                    draft = {key: overrides[question.prompt.strip('= \n')]
                             for key, question, _ in ui.fields
                             if question.prompt.strip('= \n') in overrides}
                    ui, files, report = gui.collect_stage(mkrun, answers, draft)
                    self.assertIsNone(ui.error)
                    if not ui.fields:
                        break
                    stages[ui.stage] = len(ui.fields)
                    answers.extend(ui.values)
                else:
                    self.fail('Forms did not reach preview')
                self.assertEqual(answers, expected_answers)
                self.assertEqual(files, expected_files)
                self.assertIsNotNone(report)
                self.assertEqual(stages['Run files'], 2)
                self.assertEqual(stages['Base cosmology'], 6)
                self.assertEqual(stages['AMR levels'], 2)
                if overrides.get('Zoom-in run?'):
                    self.assertEqual(stages['Zoom region'], 7)
                if overrides.get('Run mode') == 'hydro':
                    self.assertEqual(stages['Hydro solver'], 6)

    def test_stage_navigation_and_invalid_pair_are_atomic(self):
        wizard = gui.RunWizard.__new__(gui.RunWizard)
        wizard.generator, wizard.answers, wizard.history = mkrun, [], []
        wizard.draft, wizard.error = {}, mock.Mock()
        wizard.refresh = mock.Mock()
        while True:
            ui, _, _ = gui.collect_stage(mkrun, wizard.answers)
            self.assertIsNone(ui.error)
            wizard.readers = {key: mock.Mock(return_value=raw) for key, _, raw in ui.fields}
            if ui.stage == 'AMR levels':
                break
            wizard.next()
        before = list(wizard.answers)
        keys = list(wizard.readers)
        wizard.readers[keys[0]].return_value = '12'
        wizard.readers[keys[1]].return_value = '8'
        wizard.next()
        self.assertEqual(wizard.answers, before)
        wizard.error.configure.assert_called()
        wizard.readers[keys[1]].return_value = '16'
        wizard.next()
        self.assertEqual(wizard.answers, before + [12, 16])
        wizard.back()
        self.assertEqual(wizard.answers, before)
        self.assertEqual(list(wizard.draft.values()), ['12', '16'])

    def test_model_switch_drops_obsolete_branch_fields(self):
        answers = []
        while True:
            ui, _, _ = gui.collect_stage(mkrun, answers)
            if ui.stage == 'Gravity / dark-energy sector':
                break
            answers.extend(ui.values)
        choice = ui.fields[0][0]
        ui, _, _ = gui.collect_stage(mkrun, answers, {choice: 'w0wa'})
        draft = {key: raw for key, _, raw in ui.fields}
        draft[('value', 'wa')] = '0.3'
        draft[choice] = 'lcdm'
        ui, _, _ = gui.collect_stage(mkrun, answers, draft)
        self.assertIsNone(ui.error)
        self.assertEqual(ui.values, ['lcdm'])
        self.assertEqual(len(ui.fields), 1)

    def test_missing_tk_is_graceful(self):
        with mock.patch.dict(sys.modules, {'tkinter': None}), \
                contextlib.redirect_stderr(io.StringIO()) as stderr:
            self.assertEqual(gui.launch(mkrun), 2)
        self.assertIn('Tkinter is not installed', stderr.getvalue())

    def test_missing_display_is_graceful(self):
        try:
            import tkinter as tk
        except ImportError:
            self.skipTest('Tkinter not installed')
        with mock.patch.object(tk, 'Tk', side_effect=tk.TclError('no display')), \
                contextlib.redirect_stderr(io.StringIO()) as stderr:
            self.assertEqual(gui.launch(mkrun), 2)
        self.assertIn('cannot open a graphical display', stderr.getvalue())


class SaveTests(unittest.TestCase):
    def test_explicit_save_exact_bytes_and_overwrite(self):
        with tempfile.TemporaryDirectory() as directory:
            path = str(Path(directory) / 'new' / 'run.nml')
            files = {path: '&RUN_PARAMS\ncosmo=.true.\n/\n'}
            expected = gui.snapshot_targets(files)
            self.assertFalse(Path(path).parent.exists())
            self.assertEqual(gui.save_preview(files, expected), [path])
            self.assertEqual(Path(path).read_text(), files[path])
            expected = gui.snapshot_targets(files)
            files[path] += '! changed\n'
            gui.save_preview(files, expected)
            self.assertEqual(Path(path).read_text(), files[path])

    def test_target_changes_refuse_entire_batch(self):
        with tempfile.TemporaryDirectory() as directory:
            first, second = (str(Path(directory) / name) for name in ('first.nml', 'second.ini'))
            files = {first: 'new first', second: 'new second'}
            expected = gui.snapshot_targets(files)
            Path(second).write_text('concurrent work')
            with self.assertRaisesRegex(ValueError, 'destination changed'):
                gui.save_preview(files, expected)
            self.assertFalse(Path(first).exists())
            self.assertEqual(Path(second).read_text(), 'concurrent work')

    def test_symlinks_and_nonfiles_refused(self):
        with tempfile.TemporaryDirectory() as directory:
            target = Path(directory) / 'config.nml'
            target.write_text('preserve')
            link = Path(directory) / 'link.nml'
            link.symlink_to(target)
            for path in (str(link), directory):
                with self.assertRaises(ValueError):
                    gui.snapshot_targets({path: 'replace'})
            self.assertEqual(target.read_text(), 'preserve')

    def test_cancel_save_and_validation_errors_never_write(self):
        wizard = gui.RunWizard.__new__(gui.RunWizard)
        wizard.root = None
        wizard.report = {'messages': []}
        wizard.messagebox = mock.Mock()
        wizard.messagebox.askyesno.return_value = False
        with tempfile.TemporaryDirectory() as directory:
            path = str(Path(directory) / 'config.nml')
            wizard.files = {path: 'preview'}
            with mock.patch.object(gui, 'save_preview') as save:
                wizard.save()
                save.assert_not_called()
            self.assertFalse(Path(path).exists())
            wizard.report['messages'] = [mkrun.rng.ValidationMsg('ERROR', 'invalid')]
            wizard.messagebox.reset_mock()
            wizard.save()
            wizard.messagebox.askyesno.assert_not_called()


class DisplayTests(unittest.TestCase):
    def test_real_widgets_next_back_advanced_preview_and_save(self):
        try:
            import tkinter as tk
            from tkinter import ttk
        except ImportError as exc:
            self.skipTest('Tkinter unavailable: {}'.format(exc))
        try:
            root = tk.Tk()
        except tk.TclError as exc:
            # The model/save suite still runs on machines without Tk/display.
            self.skipTest('Graphical display unavailable: {}'.format(exc))
        root.withdraw()
        try:
            dialogs = (mock.Mock(), mock.Mock())
            wizard = gui.RunWizard(root, mkrun, tk, ttk, dialogs)
            with tempfile.TemporaryDirectory() as directory:
                name_key = ('value', 'Run name (used as file/dir prefix)')
                wizard.controls[name_key].set('gui_smoke')
                wizard.controls[('value', 'Output directory')].set(directory)
                wizard.next()
                wizard.back()
                self.assertEqual(wizard.controls[name_key].get(), 'gui_smoke')
                wizard.next()
                for _ in range(100):
                    root.update_idletasks()
                    if wizard.report is not None:
                        break
                    for key, question, _ in wizard.stage_ui.fields:
                        if question.kind == 'choice' and 'Run mode' in question.prompt:
                            wizard.controls[key].set('hydro')
                        if question.kind == 'bool' and ('full parameter editor' in question.prompt
                                                        or 'Zoom-in run' in question.prompt):
                            wizard.controls[key].set(True)
                    wizard.update_fields()
                    if wizard.stage_ui.stage == 'Base cosmology':
                        self.assertEqual(len(wizard.controls), 6)
                    if wizard.stage_ui.stage == 'Zoom region':
                        self.assertEqual(len(wizard.controls), 7)
                    wizard.next()
                self.assertIsNotNone(wizard.report)
                self.assertFalse(list(Path(directory).iterdir()))
                self.assertTrue(wizard.report['values']['hydro'])
                self.assertEqual(str(wizard.save_button['state']), 'normal')

                def widgets(parent):
                    for child in parent.winfo_children():
                        yield child
                        yield from widgets(child)

                preview_widgets = [child for child in widgets(wizard.frame) if isinstance(child, tk.Text)]
                self.assertEqual(len(preview_widgets), len(wizard.files) + 1)
                self.assertTrue(all(child['state'] == 'disabled' for child in preview_widgets))
                dialogs[1].askyesno.return_value = False
                wizard.save()
                self.assertFalse(list(Path(directory).iterdir()))
                dialogs[1].askyesno.return_value = True
                wizard.save()
                self.assertEqual({path: Path(path).read_text() for path in wizard.files}, wizard.files)
                dialogs[1].showerror.assert_not_called()
        finally:
            root.destroy()


if __name__ == '__main__':
    unittest.main()
