# -*- coding: utf-8 -*-
"""
3D GRID SPLITTER - GUI
======================
grid_splitter.py 의 기능을 창에서 사용할 수 있게 만든 Tkinter 인터페이스.
PyInstaller 로 EXE 를 만들면 파이썬 설치 없이 실행할 수 있다 (build_exe.bat 참고).
"""

import os
import sys
import queue
import threading
import traceback
import subprocess

import tkinter as tk
from tkinter import ttk, filedialog, messagebox

# PyInstaller onefile 로 묶였을 때도 모듈을 찾도록
if getattr(sys, "frozen", False):
    sys.path.insert(0, os.path.dirname(sys.executable))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import grid_splitter as gs

APP_TITLE = "3D GRID SPLITTER"
PAD = 8


class Cancelled(Exception):
    pass


class App(ttk.Frame):
    def __init__(self, master):
        super().__init__(master, padding=PAD)
        self.grid(sticky="nsew")
        master.columnconfigure(0, weight=1)
        master.rowconfigure(0, weight=1)
        self.columnconfigure(0, weight=1)

        self.q = queue.Queue()
        self.worker = None
        self.cancel = threading.Event()
        self.bbox = None            # 선택한 모델의 (min, max)
        self.last_out = None

        self._build()
        self.after(80, self._drain)

    # ----------------------------- 화면 구성 -----------------------------
    def _build(self):
        r = 0
        # --- 입력 / 출력 ---
        box = ttk.LabelFrame(self, text=" 파일 ", padding=PAD)
        box.grid(row=r, column=0, sticky="ew", pady=(0, PAD)); r += 1
        box.columnconfigure(1, weight=1)

        ttk.Label(box, text="입력 모델").grid(row=0, column=0, sticky="w")
        self.v_in = tk.StringVar()
        ttk.Entry(box, textvariable=self.v_in).grid(row=0, column=1, sticky="ew", padx=6)
        ttk.Button(box, text="찾아보기…", command=self._pick_in).grid(row=0, column=2)

        ttk.Label(box, text="저장 폴더").grid(row=1, column=0, sticky="w", pady=(6, 0))
        self.v_out = tk.StringVar()
        ttk.Entry(box, textvariable=self.v_out).grid(row=1, column=1, sticky="ew",
                                                     padx=6, pady=(6, 0))
        ttk.Button(box, text="찾아보기…", command=self._pick_out).grid(row=1, column=2,
                                                                    pady=(6, 0))
        self.v_info = tk.StringVar(value="모델을 선택하면 크기와 예상 조각 수가 표시됩니다.")
        ttk.Label(box, textvariable=self.v_info, foreground="#41546b").grid(
            row=2, column=0, columnspan=3, sticky="w", pady=(8, 0))

        # --- 절단 크기 ---
        box = ttk.LabelFrame(self, text=" 절단 크기 (mm) ", padding=PAD)
        box.grid(row=r, column=0, sticky="ew", pady=(0, PAD)); r += 1
        self.v_x = tk.StringVar(value="200")
        self.v_y = tk.StringVar(value="200")
        self.v_z = tk.StringVar(value="200")
        for i, (lab, var) in enumerate([("X", self.v_x), ("Y", self.v_y), ("Z", self.v_z)]):
            ttk.Label(box, text=lab).grid(row=0, column=i * 2, sticky="e", padx=(12 if i else 0, 4))
            e = ttk.Entry(box, textvariable=var, width=9, justify="right")
            e.grid(row=0, column=i * 2 + 1, sticky="w")
            var.trace_add("write", lambda *a: self._refresh_info())
        ttk.Button(box, text="X값을 Y·Z에 복사", command=self._same).grid(
            row=0, column=6, padx=(16, 0))

        # --- 옵션 ---
        box = ttk.LabelFrame(self, text=" 옵션 ", padding=PAD)
        box.grid(row=r, column=0, sticky="ew", pady=(0, PAD)); r += 1
        self.v_fmt = tk.StringVar(value="stl")
        ttk.Label(box, text="저장 형식").grid(row=0, column=0, sticky="w")
        ttk.Combobox(box, textvariable=self.v_fmt, values=["stl", "obj"], width=6,
                     state="readonly").grid(row=0, column=1, sticky="w", padx=(6, 18))

        self.v_repair = tk.BooleanVar(value=True)
        self.v_norm = tk.BooleanVar(value=True)
        self.v_layer = tk.BooleanVar(value=True)
        self.v_expl = tk.BooleanVar(value=True)
        self.v_csv = tk.BooleanVar(value=True)
        self.v_zip = tk.BooleanVar(value=False)
        ttk.Checkbutton(box, text="원본 구멍 메우기", variable=self.v_repair).grid(
            row=0, column=2, sticky="w")
        ttk.Checkbutton(box, text="법선 방향 자동 정리", variable=self.v_norm).grid(
            row=0, column=3, sticky="w", padx=(18, 0))
        ttk.Checkbutton(box, text="층별 평면도 PNG", variable=self.v_layer).grid(
            row=1, column=0, columnspan=2, sticky="w", pady=(6, 0))
        ttk.Checkbutton(box, text="3D 분해 조립도 PNG", variable=self.v_expl).grid(
            row=1, column=2, sticky="w", pady=(6, 0))
        ttk.Checkbutton(box, text="조각 목록 index.csv", variable=self.v_csv).grid(
            row=2, column=0, columnspan=2, sticky="w", pady=(6, 0))
        ttk.Checkbutton(box, text="결과를 ZIP 으로 묶기", variable=self.v_zip).grid(
            row=2, column=2, sticky="w", pady=(6, 0))

        ttk.Label(box, text="분해도 벌림").grid(row=3, column=0, sticky="w", pady=(8, 0))
        self.v_explode = tk.DoubleVar(value=0.25)
        ttk.Scale(box, from_=0.0, to=1.0, variable=self.v_explode,
                  orient="horizontal", length=180).grid(row=3, column=1, columnspan=2,
                                                        sticky="w", padx=6, pady=(8, 0))

        # --- 실행 ---
        bar = ttk.Frame(self)
        bar.grid(row=r, column=0, sticky="ew", pady=(0, PAD)); r += 1
        bar.columnconfigure(2, weight=1)
        self.btn_run = ttk.Button(bar, text="자르기 시작", command=self._start)
        self.btn_run.grid(row=0, column=0)
        self.btn_stop = ttk.Button(bar, text="중지", command=self._stop, state="disabled")
        self.btn_stop.grid(row=0, column=1, padx=6)
        self.pb = ttk.Progressbar(bar, mode="determinate", maximum=1000)
        self.pb.grid(row=0, column=2, sticky="ew", padx=6)
        self.btn_open = ttk.Button(bar, text="결과 폴더 열기", command=self._open_out,
                                   state="disabled")
        self.btn_open.grid(row=0, column=3)

        # --- 로그 ---
        box = ttk.LabelFrame(self, text=" 진행 상황 ", padding=4)
        box.grid(row=r, column=0, sticky="nsew"); r += 1
        self.rowconfigure(r - 1, weight=1)
        box.columnconfigure(0, weight=1); box.rowconfigure(0, weight=1)
        self.txt = tk.Text(box, height=16, wrap="none", font=("Consolas", 9),
                           background="#12181f", foreground="#dfe7ef",
                           insertbackground="#dfe7ef", relief="flat")
        self.txt.grid(row=0, column=0, sticky="nsew")
        sb = ttk.Scrollbar(box, command=self.txt.yview)
        sb.grid(row=0, column=1, sticky="ns")
        self.txt.configure(yscrollcommand=sb.set, state="disabled")

        self._log("3D GRID SPLITTER — 모델을 격자로 잘라 1.stl, 2.stl … 로 저장하고")
        self._log("각 조각의 위치를 보여주는 층별 지도를 만듭니다.")
        self._log("번호 순서: 아래층부터 → 도면 위쪽 행부터 → 왼쪽 열부터")
        self._log("")

    # ----------------------------- 동작 -----------------------------
    def _pick_in(self):
        p = filedialog.askopenfilename(
            title="3D 모델 선택",
            filetypes=[("3D 모델", "*.stl *.obj"), ("STL", "*.stl"),
                       ("OBJ", "*.obj"), ("모든 파일", "*.*")])
        if not p:
            return
        self.v_in.set(p)
        if not self.v_out.get():
            self.v_out.set(os.path.join(os.path.dirname(p),
                                        os.path.splitext(os.path.basename(p))[0] + "_pieces"))
        self.bbox = None
        self.v_info.set("모델 크기를 읽는 중…")
        threading.Thread(target=self._probe, args=(p,), daemon=True).start()

    def _pick_out(self):
        p = filedialog.askdirectory(title="저장 폴더 선택")
        if p:
            self.v_out.set(p)

    def _same(self):
        self.v_y.set(self.v_x.get()); self.v_z.set(self.v_x.get())

    def _probe(self, path):
        try:
            V, F = gs.load_mesh(path)
            self.q.put(("bbox", (V.min(0), V.max(0), len(V), len(F))))
        except Exception as ex:
            self.q.put(("info", f"파일을 읽지 못했습니다: {ex}"))

    def _sizes(self):
        try:
            sx, sy, sz = (float(self.v_x.get()), float(self.v_y.get()), float(self.v_z.get()))
        except ValueError:
            return None
        if min(sx, sy, sz) <= 0:
            return None
        return sx, sy, sz

    def _refresh_info(self):
        if self.bbox is None:
            return
        mn, mx, nv, nf = self.bbox
        d = mx - mn
        s = self._sizes()
        if s is None:
            self.v_info.set(f"크기 {d[0]:.1f} × {d[1]:.1f} × {d[2]:.1f} mm  |  "
                            f"절단 크기를 숫자로 입력하세요")
            return
        n = [max(1, int(np.ceil(d[i] / s[i]))) for i in range(3)]
        self.v_info.set(
            f"크기 {d[0]:.1f} × {d[1]:.1f} × {d[2]:.1f} mm  |  정점 {nv:,} / 면 {nf:,}  |  "
            f"격자 {n[0]}열 × {n[1]}행 × {n[2]}층 → 최대 {n[0]*n[1]*n[2]}조각")

    def _log(self, msg=""):
        self.txt.configure(state="normal")
        self.txt.insert("end", str(msg) + "\n")
        self.txt.see("end")
        self.txt.configure(state="disabled")

    def _start(self):
        if self.worker and self.worker.is_alive():
            return
        src = self.v_in.get().strip()
        out = self.v_out.get().strip()
        s = self._sizes()
        if not src or not os.path.exists(src):
            messagebox.showwarning(APP_TITLE, "입력 모델 파일을 선택하세요."); return
        if not out:
            messagebox.showwarning(APP_TITLE, "저장 폴더를 선택하세요."); return
        if s is None:
            messagebox.showwarning(APP_TITLE, "절단 크기는 0보다 큰 숫자여야 합니다."); return

        self.cancel.clear()
        self.btn_run.configure(state="disabled")
        self.btn_stop.configure(state="normal")
        self.btn_open.configure(state="disabled")
        self.pb["value"] = 0
        self._log("─" * 70)

        opts = dict(cut_x=s[0], cut_y=s[1], cut_z=s[2],
                    output_format=self.v_fmt.get(),
                    repair_source_holes=self.v_repair.get(),
                    fix_normals=self.v_norm.get(),
                    layer_maps=self.v_layer.get(),
                    exploded_map=self.v_expl.get(),
                    index_csv=self.v_csv.get(),
                    make_zip=self.v_zip.get(),
                    explode=float(self.v_explode.get()))
        self.worker = threading.Thread(target=self._run, args=(src, out, opts), daemon=True)
        self.worker.start()

    def _stop(self):
        self.cancel.set()
        self._log("중지 요청… 현재 단계가 끝나면 멈춥니다.")

    def _run(self, src, out, opts):
        def log(m=""):
            self.q.put(("log", m))

        def prog(f):
            if self.cancel.is_set():
                raise Cancelled()
            self.q.put(("prog", float(f)))

        try:
            res = gs.run_split(src, out, log=log, progress=prog, **opts)
            self.q.put(("done", (out, res)))
        except Cancelled:
            self.q.put(("log", "사용자가 중지했습니다."))
            self.q.put(("fail", None))
        except Exception as ex:
            self.q.put(("log", "오류: " + str(ex)))
            self.q.put(("log", traceback.format_exc()))
            self.q.put(("fail", str(ex)))

    def _drain(self):
        try:
            while True:
                kind, payload = self.q.get_nowait()
                if kind == "log":
                    self._log(payload)
                elif kind == "prog":
                    self.pb["value"] = max(0, min(1000, int(payload * 1000)))
                elif kind == "bbox":
                    self.bbox = payload
                    self._refresh_info()
                elif kind == "info":
                    self.v_info.set(payload)
                elif kind == "done":
                    out, res = payload
                    self.last_out = out
                    self.btn_run.configure(state="normal")
                    self.btn_stop.configure(state="disabled")
                    self.btn_open.configure(state="normal")
                    self.pb["value"] = 1000
                    messagebox.showinfo(
                        APP_TITLE,
                        f"완료!\n\n조각 {len(res['pieces'])}개 / 층 {res['n_layers']}개\n"
                        f"지도 {len(res['maps'])}장\n\n저장 위치:\n{out}")
                elif kind == "fail":
                    self.btn_run.configure(state="normal")
                    self.btn_stop.configure(state="disabled")
                    if payload:
                        messagebox.showerror(APP_TITLE, f"작업을 마치지 못했습니다.\n\n{payload}")
        except queue.Empty:
            pass
        self.after(80, self._drain)

    def _open_out(self):
        p = self.last_out or self.v_out.get()
        if not p or not os.path.isdir(p):
            return
        try:
            if os.name == "nt":
                os.startfile(p)                       # noqa: S606
            elif sys.platform == "darwin":
                subprocess.Popen(["open", p])
            else:
                subprocess.Popen(["xdg-open", p])
        except Exception as ex:
            messagebox.showinfo(APP_TITLE, f"폴더를 열지 못했습니다: {ex}")


def main():
    if os.name == "nt":                               # 고해상도 모니터에서 또렷하게
        try:
            import ctypes
            ctypes.windll.shcore.SetProcessDpiAwareness(1)
        except Exception:
            pass
    root = tk.Tk()
    root.title(APP_TITLE)
    root.geometry("860x740")
    root.minsize(720, 600)
    try:
        ttk.Style().theme_use("vista" if os.name == "nt" else "clam")
    except Exception:
        pass
    App(root)
    root.mainloop()


if __name__ == "__main__":
    main()
