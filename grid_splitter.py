# -*- coding: utf-8 -*-
"""
3D GRID SPLITTER - core engine
==============================
큰 3D 모델(.stl/.obj)을 격자로 잘라 방수(watertight) 조각으로 저장한다.

원본(grid cut.py) 대비 추가된 기능
  1) 파일명을 1.stl, 2.stl, 3.stl ... 처럼 숫자 통번호로 저장
     번호 순서: 아래층(Z 작은 쪽)부터 -> 각 층에서 위쪽 행(Y 큰 쪽)부터 -> 왼쪽 열(X 작은 쪽)부터
  2) 조각 위치 지도 생성
     - layer_01.png ... : 각 층을 위에서 내려다본 평면도 (실제 조각 단면 + 번호)
     - layers_overview.png : 모든 층을 한 장에 모은 요약 도면
     - exploded_3d.png : 3D 분해 조립도
     - index.csv : 번호 / 층 / 행 / 열 / 좌표범위 / 면 개수 목록
  3) GUI(gui.py) 및 CLI에서 모두 호출할 수 있도록 함수화 (진행률 콜백 지원)

의존성: numpy (필수), matplotlib (지도 생성 시), scipy (있으면 사용, 없으면 자체 구현으로 대체)
"""

import os
import csv
import glob
import math
import struct
import zipfile
from collections import defaultdict

import numpy as np

# scipy 가 있으면 cKDTree 를 쓰고, 없으면 numpy 격자 해싱으로 대체한다.
try:
    from scipy.spatial import cKDTree as _cKDTree
    _HAS_SCIPY = True
except Exception:          # pragma: no cover
    _cKDTree = None
    _HAS_SCIPY = False


def _pairs_within_numpy(V, r):
    """scipy 없이 반지름 r 안의 점 쌍을 찾는다 (격자 해싱). 반환: (M,2) 인덱스 배열."""
    if len(V) < 2:
        return np.zeros((0, 2), np.int64)
    c = np.floor(V / r).astype(np.int64)
    c = c - c.min(0) + 1                     # 이웃(-1) 조회해도 음수가 되지 않도록 +1
    if c.max() + 1 >= (1 << 20):             # 격자 범위 초과 -> 쌍 없음으로 처리
        return np.zeros((0, 2), np.int64)
    def keyof(a):
        return (a[:, 0] << 40) | (a[:, 1] << 20) | a[:, 2]
    key = keyof(c)
    uk, inv = np.unique(key, return_inverse=True)
    order = np.argsort(inv, kind="stable")
    counts = np.bincount(inv, minlength=len(uk))
    starts = np.concatenate([[0], np.cumsum(counts)[:-1]])
    offs = [(dx, dy, dz)
            for dx in (0, 1) for dy in (-1, 0, 1) for dz in (-1, 0, 1)
            if (dx, dy, dz) >= (0, 0, 0)]    # 쌍은 대칭이므로 절반만 조회
    out = []
    for off in offs:
        qk = keyof(c + np.array(off, np.int64))
        pos = np.searchsorted(uk, qk)
        np.clip(pos, 0, len(uk) - 1, out=pos)
        hit = uk[pos] == qk
        if not hit.any():
            continue
        pi = np.nonzero(hit)[0]
        gi = pos[hit]
        cnt = counts[gi]
        total = int(cnt.sum())
        if total == 0:
            continue
        rep_i = np.repeat(pi, cnt)
        ramp = np.arange(total) - np.repeat(np.cumsum(cnt) - cnt, cnt)
        rep_j = order[np.repeat(starts[gi], cnt) + ramp]
        m = rep_i < rep_j
        if m.any():
            out.append(np.column_stack([rep_i[m], rep_j[m]]))
    if not out:
        return np.zeros((0, 2), np.int64)
    P = np.vstack(out)
    d = np.linalg.norm(V[P[:, 0]] - V[P[:, 1]], axis=1)
    P = P[d <= r]
    if len(P) == 0:
        return np.zeros((0, 2), np.int64)
    return np.unique(P, axis=0)


def _query_pairs(V, r):
    if _HAS_SCIPY:
        return _cKDTree(V).query_pairs(r=r, output_type="ndarray")
    return _pairs_within_numpy(V, r)


# ============================ 2D 삼각분할 (구멍 포함) ============================
def _signed_area(P):
    x = P[:, 0]; y = P[:, 1]
    return 0.5 * np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y)

def _pip(p, poly):
    x, y = p; inside = False; m = len(poly); j = m - 1
    for i in range(m):
        xi, yi = poly[i]; xj, yj = poly[j]
        if ((yi > y) != (yj > y)) and (x < (xj - xi) * (y - yi) / (yj - yi + 1e-30) + xi):
            inside = not inside
        j = i
    return inside

def _point_in_tri(p, a, b, c):
    d1 = (b[0]-a[0])*(p[1]-a[1]) - (b[1]-a[1])*(p[0]-a[0])
    d2 = (c[0]-b[0])*(p[1]-b[1]) - (c[1]-b[1])*(p[0]-b[0])
    d3 = (a[0]-c[0])*(p[1]-c[1]) - (a[1]-c[1])*(p[0]-c[0])
    neg = (d1 < 0) or (d2 < 0) or (d3 < 0)
    pos = (d1 > 0) or (d2 > 0) or (d3 > 0)
    return not (neg and pos)

def _ear_clip(poly_idx, P2):
    idx = list(poly_idx)
    if len(idx) < 3: return []
    if _signed_area(P2[idx]) < 0: idx = idx[::-1]
    pts = P2; tris = []
    def cross(a, b, c):
        return (pts[b,0]-pts[a,0])*(pts[c,1]-pts[a,1]) - (pts[b,1]-pts[a,1])*(pts[c,0]-pts[a,0])
    guard = 0
    while len(idx) > 3 and guard < 5*len(poly_idx) + 10000:
        guard += 1; m = len(idx); found = False
        for i in range(m):
            a, b, c = idx[(i-1)%m], idx[i], idx[(i+1)%m]
            if cross(a, b, c) <= 1e-14: continue
            pa, pb, pc = pts[a], pts[b], pts[c]; ok = True
            for j in idx:
                if j in (a, b, c): continue
                if _point_in_tri(pts[j], pa, pb, pc): ok = False; break
            if ok: tris.append((a, b, c)); del idx[i]; found = True; break
        if not found:                      # 안전장치: 무한루프 방지
            best = None; bestv = -1e30
            for i in range(m):
                a, b, c = idx[(i-1)%m], idx[i], idx[(i+1)%m]; cv = cross(a, b, c)
                if cv > bestv: bestv = cv; best = i
            if best is None: break
            i = best; a, b, c = idx[(i-1)%m], idx[i], idx[(i+1)%m]
            tris.append((a, b, c)); del idx[i]
    if len(idx) == 3: tris.append((idx[0], idx[1], idx[2]))
    return tris

def _bridge(outer, hole, P2):
    hi = max(range(len(hole)), key=lambda i: P2[hole[i], 0]); M = P2[hole[hi]]
    m = len(outer); bestx = np.inf; edge = None; Ipt = None
    for i in range(m):
        a = P2[outer[i]]; b = P2[outer[(i+1)%m]]; ya, yb = a[1], b[1]
        if (ya > M[1]) == (yb > M[1]): continue
        t = (M[1]-ya)/(yb-ya); xi = a[0] + t*(b[0]-a[0])
        if xi >= M[0]-1e-12 and xi < bestx:
            bestx = xi; edge = (i, (i+1)%m); Ipt = np.array([xi, M[1]])
    if edge is None:
        vis = min(range(m), key=lambda i: np.sum((P2[outer[i]]-M)**2))
    else:
        P = edge[0] if P2[outer[edge[0]],0] > P2[outer[edge[1]],0] else edge[1]
        Pp = P2[outer[P]]; cand = P; bestang = -np.inf
        for vi in range(m):
            v = outer[vi]
            if vi == P: continue
            pv = P2[v]
            if _point_in_tri(pv, M, Ipt, Pp):
                d = pv - M; ang = d[0]/(np.linalg.norm(d)+1e-30)
                if ang > bestang: bestang = ang; cand = vi
        vis = cand
    hrot = hole[hi:] + hole[:hi]
    return outer[:vis+1] + hrot + [hrot[0]] + outer[vis:]

def triangulate_with_holes(loops, P2):
    """loops: P2 인덱스 리스트들. even-odd 중첩으로 외곽/구멍 자동 판별 후 삼각분할."""
    if not loops: return []
    areas = [_signed_area(P2[l]) for l in loops]
    reps  = [P2[l[0]] for l in loops]
    depth = [0]*len(loops)
    for i in range(len(loops)):
        for j in range(len(loops)):
            if i == j: continue
            if _pip(reps[i], P2[loops[j]]): depth[i] += 1
    L = [list(l) for l in loops]
    outers = [i for i in range(len(loops)) if depth[i] % 2 == 0]
    holes  = [i for i in range(len(loops)) if depth[i] % 2 == 1]
    for i in outers:
        if areas[i] < 0: L[i] = L[i][::-1]   # CCW
    for i in holes:
        if areas[i] > 0: L[i] = L[i][::-1]   # CW
    tris = []
    for oi in outers:
        opoly = P2[L[oi]]
        mine = [hi for hi in holes if _pip(reps[hi], opoly) and depth[hi] == depth[oi]+1]
        merged = list(L[oi])
        for hi in sorted(mine, key=lambda h: -max(P2[L[h]][:, 0])):
            merged = _bridge(merged, L[hi], P2)
        tris += _ear_clip(merged, P2)
    return tris

# ================================ 메시 IO ================================
def load_mesh(path):
    ext = os.path.splitext(path)[1].lower()
    if ext == ".stl": return _load_stl(path)
    if ext == ".obj": return _load_obj(path)
    raise ValueError("지원하지 않는 형식: " + ext)

def save_mesh(path, V, F):
    ext = os.path.splitext(path)[1].lower()
    if ext == ".stl": return _save_stl_binary(path, V, F)
    if ext == ".obj": return _save_obj(path, V, F)
    raise ValueError("지원하지 않는 형식: " + ext)

def _load_stl(path):
    with open(path, "rb") as f:
        head = f.read(5); f.seek(0); data = f.read()
    if head == b"solid" and b"facet" in data[:2048]:
        return _load_stl_ascii(data.decode("utf-8", "ignore"))
    return _load_stl_binary(data)

def _load_stl_binary(data):
    n = struct.unpack("<I", data[80:84])[0]
    rec = np.frombuffer(data[84:84+50*n], dtype=np.uint8).reshape(n, 50)
    floats = rec[:, :48].copy().view("<f4").reshape(n, 12)
    tris = floats[:, 3:12].reshape(n, 3, 3).astype(np.float64)
    V = tris.reshape(-1, 3); F = np.arange(len(V)).reshape(-1, 3)
    return weld(V, F)

def _load_stl_ascii(text):
    verts = []
    for line in text.splitlines():
        s = line.strip()
        if s.startswith("vertex"): verts.append([float(x) for x in s.split()[1:4]])
    V = np.array(verts, float); F = np.arange(len(V)).reshape(-1, 3)
    return weld(V, F)

def _load_obj(path):
    verts, faces = [], []
    with open(path) as f:
        for line in f:
            if line.startswith("v "):
                verts.append([float(x) for x in line.split()[1:4]])
            elif line.startswith("f "):
                idx = [int(p.split("/")[0]) - 1 for p in line.split()[1:]]
                for k in range(1, len(idx)-1): faces.append([idx[0], idx[k], idx[k+1]])
    return np.array(verts, float), np.array(faces, np.int64)

def _save_stl_binary(path, V, F):
    tris = V[F]; e1 = tris[:,1]-tris[:,0]; e2 = tris[:,2]-tris[:,0]
    nor = np.cross(e1, e2); ln = np.linalg.norm(nor, axis=1, keepdims=True)
    ln[ln == 0] = 1.0; nor /= ln
    with open(path, "wb") as f:
        f.write(b"\0"*80); f.write(struct.pack("<I", len(F)))
        out = bytearray()
        for i in range(len(F)):
            out += struct.pack("<3f", *nor[i])
            for v in tris[i]: out += struct.pack("<3f", *v)
            out += struct.pack("<H", 0)
        f.write(out)

def _save_obj(path, V, F):
    with open(path, "w") as f:
        for v in V: f.write(f"v {v[0]} {v[1]} {v[2]}\n")
        for t in F: f.write(f"f {t[0]+1} {t[1]+1} {t[2]+1}\n")

# ================================ 메시 유틸 ================================
def _unique_rows(a):
    order = np.lexsort(a.T[::-1]); a_s = a[order]
    diff = np.ones(len(a), bool); diff[1:] = np.any(a_s[1:] != a_s[:-1], axis=1)
    uid = np.cumsum(diff) - 1; inv = np.empty(len(a), np.int64); inv[order] = uid
    return None, inv, order[diff]

def weld(V, F, tol=1e-7):
    if len(V) == 0: return V, F
    diag = np.linalg.norm(V.max(0) - V.min(0)) or 1.0
    q = np.round(V / (tol * diag)).astype(np.int64)
    _, inv, idx = _unique_rows(q); Vn = V[idx]; Fn = inv[F]
    good = (Fn[:,0] != Fn[:,1]) & (Fn[:,1] != Fn[:,2]) & (Fn[:,0] != Fn[:,2])
    return Vn, Fn[good]

def robust_weld(V, F, tol):
    if len(V) == 0: return V, F
    pairs = _query_pairs(V, tol)
    parent = np.arange(len(V))
    def find(x):
        r = x
        while parent[r] != r: r = parent[r]
        while parent[x] != r: parent[x], x = r, parent[x]
        return r
    for i, j in pairs:
        ri, rj = find(i), find(j)
        if ri != rj: parent[max(ri, rj)] = min(ri, rj)
    root = np.array([find(i) for i in range(len(V))])
    uniq, inv = np.unique(root, return_inverse=True); Vn = V[uniq]; Fn = inv[F]
    good = (Fn[:,0] != Fn[:,1]) & (Fn[:,1] != Fn[:,2]) & (Fn[:,0] != Fn[:,2])
    return Vn, Fn[good]

def watertight_report(V, F):
    if len(F) == 0:
        return {"watertight": False, "faces": 0, "boundary_edges": 0, "nonmanifold": 0}
    e = np.sort(np.vstack([F[:,[0,1]], F[:,[1,2]], F[:,[2,0]]]), axis=1)
    _, inv, _ = _unique_rows(e); counts = np.bincount(inv)
    b = int(np.sum(counts == 1)); nm = int(np.sum(counts > 2))
    return {"watertight": b == 0 and nm == 0, "faces": len(F),
            "boundary_edges": b, "nonmanifold": nm}

def _project(Vp, n, p0):
    n = np.asarray(n, float); n = n / (np.linalg.norm(n) or 1.0)
    u = np.cross(n, [1., 0, 0])
    if np.linalg.norm(u) < 1e-6: u = np.cross(n, [0, 1., 0])
    u /= np.linalg.norm(u); w = np.cross(n, u)
    return np.column_stack([(Vp - p0) @ u, (Vp - p0) @ w])

# ================== 단면 세그먼트 → 닫힌 루프 (분기점 견고 + gap 보정) ==================
def _chain_loops_from_segs(segs, n, p0, diag):
    """방향성 절단 세그먼트를 닫힌 루프로 조립.
       분기점(한 점에서 여러 갈래)에서 턴 각도가 가장 일관된 방향을 따라가
       외곽을 정확히 추적한다(스퍼/자기접촉에 강함). 작은 틈은 시작점으로 닫는다."""
    pts = segs.reshape(-1, 3); scale = 1e-6 * diag
    q = np.round(pts / scale).astype(np.int64)
    uq, inv = np.unique(q, axis=0, return_inverse=True)
    K = len(uq); Vp = np.zeros((K, 3)); cnt = np.zeros(K)
    np.add.at(Vp, inv, pts); np.add.at(cnt, inv, 1.0); Vp /= cnt[:, None]
    e = inv.reshape(-1, 2)
    P2 = _project(Vp, n, p0)
    loops = _assemble_loops_geometric(e, P2, gap_tol=0.02 * diag)
    return loops, Vp

def _assemble_loops_geometric(edges, P2, gap_tol):
    from collections import defaultdict
    adj = defaultdict(list)
    for a, b in edges:
        if a != b: adj[a].append(b); adj[b].append(a)
    used = set()
    def ekey(u, v): return (u, v) if u < v else (v, u)
    loops = []
    starts = sorted(adj.keys(), key=lambda v: (len(adj[v]) % 2 == 0, len(adj[v])))
    ne = len(edges)
    for s in starts:
        for first in adj[s]:
            if ekey(s, first) in used: continue
            loop = [s]; prev = s; cur = first; used.add(ekey(s, first))
            while cur != s:
                loop.append(cur)
                inc = P2[cur] - P2[prev]; ia = np.arctan2(inc[1], inc[0])
                best = None; bestturn = None
                for w in adj[cur]:
                    if ekey(cur, w) in used: continue
                    d = P2[w] - P2[cur]; a = np.arctan2(d[1], d[0])
                    turn = (a - ia) % (2 * np.pi)      # 가장 작은 CCW 턴 = 외곽 우선
                    if bestturn is None or turn < bestturn: bestturn = turn; best = w
                if best is None: break                 # dead-end: 시작점으로 암묵적 닫힘
                used.add(ekey(cur, best)); prev = cur; cur = best
                if len(loop) > ne + 5: break
            if len(loop) >= 3: loops.append(loop)
    return loops

# ===================== 반평면 절단 + 즉시 캡 =====================
def clip_and_cap(V, F, n, p0, tol, diag):
    """(x-p0)·n >= 0 인 쪽만 남기고, 그 절단으로 생긴 단면을 즉시 막는다."""
    if len(F) == 0: return V, F
    n = np.asarray(n, float); n = n / (np.linalg.norm(n) or 1.0); p0 = np.asarray(p0, float)
    d = (V - p0) @ n; d[np.abs(d) < tol] = 0.0
    tri = V[F]; dd = d[F]; inside = dd >= 0; cnt = inside.sum(1)
    parts = []; segs = []
    if (cnt == 3).any(): parts.append(tri[cnt == 3])
    def interp(A, B, dA, dB):
        t = (dA / (dA - dB))[:, None]; return A + t * (B - A)
    one = cnt == 1
    if one.any():
        T = tri[one]; D = dd[one]; IN = inside[one]
        roll = np.argmax(IN, 1); idx = (np.arange(3)[None] + roll[:, None]) % 3
        T = np.take_along_axis(T, idx[:, :, None].repeat(3, 2), 1); D = np.take_along_axis(D, idx, 1)
        A, B, C = T[:,0], T[:,1], T[:,2]; dA, dB, dC = D[:,0], D[:,1], D[:,2]
        P = interp(A, B, dA, dB); Q = interp(A, C, dA, dC)
        parts.append(np.stack([A, P, Q], 1)); segs.append(np.stack([Q, P], 1))
    two = cnt == 2
    if two.any():
        T = tri[two]; D = dd[two]; IN = inside[two]
        roll = np.argmin(IN, 1); idx = (np.arange(3)[None] + roll[:, None]) % 3
        T = np.take_along_axis(T, idx[:, :, None].repeat(3, 2), 1); D = np.take_along_axis(D, idx, 1)
        A, B, C = T[:,0], T[:,1], T[:,2]; dA, dB, dC = D[:,0], D[:,1], D[:,2]
        P = interp(A, B, dA, dB); Q = interp(A, C, dA, dC)
        parts.append(np.stack([B, C, Q], 1)); parts.append(np.stack([B, Q, P], 1))
        segs.append(np.stack([P, Q], 1))
    if not parts: return np.zeros((0, 3)), np.zeros((0, 3), int)
    wall = np.concatenate(parts, 0)
    cap_tris = []
    if segs:
        S = np.concatenate(segs, 0)
        loops, Vp = _chain_loops_from_segs(S, n, p0, diag)
        if loops:
            P2 = _project(Vp, n, p0)
            tl = triangulate_with_holes(loops, P2)
            if tl:
                Tl = np.array(tl, int); capV = Vp[Tl]
                nrm = np.cross(capV[:,1]-capV[:,0], capV[:,2]-capV[:,0])
                flip = (nrm @ (-n)) < 0
                capV[flip] = capV[flip][:, ::-1]
                cap_tris.append(capV)
    allt = wall if not cap_tris else np.concatenate([wall] + cap_tris, 0)
    Vc = allt.reshape(-1, 3); Fc = np.arange(len(Vc)).reshape(-1, 3)
    Vc, Fc = robust_weld(Vc, Fc, tol=5e-6 * diag)
    return Vc, Fc

def slab(V, F, axis, lo, hi, tol, diag):
    nplus = np.zeros(3); nplus[axis] = 1.0
    plo = np.zeros(3); plo[axis] = lo
    V, F = clip_and_cap(V, F, nplus, plo, tol, diag)
    if len(F) == 0: return V, F
    phi = np.zeros(3); phi[axis] = hi
    V, F = clip_and_cap(V, F, -nplus, phi, tol, diag)
    return V, F

# =============================== 구멍 메우기(선택) ===============================
def fill_holes(V, F, max_iter=8):
    """원본 스캔의 기존 경계 구멍을 루프별 최적평면 삼각분할로 메운다."""
    if len(F) == 0: return V, F
    diag = np.linalg.norm(V.max(0) - V.min(0)) or 1.0
    for _ in range(max_iter):
        he = np.vstack([F[:,[0,1]], F[:,[1,2]], F[:,[2,0]]]); key = np.sort(he, axis=1)
        _, inv, _ = _unique_rows(key); counts = np.bincount(inv)
        bnd = he[counts[inv] == 1]
        if len(bnd) == 0: break
        succ = defaultdict(list)
        for a, b in bnd:
            if a != b: succ[a].append(b)
        ptr = defaultdict(int); loops = []
        for s in list(succ.keys()):
            while ptr[s] < len(succ[s]):
                loop = [s]; cur = s; ok = True
                while True:
                    if ptr[cur] >= len(succ[cur]): ok = False; break
                    nxt = succ[cur][ptr[cur]]; ptr[cur] += 1
                    if nxt == s: break
                    loop.append(nxt); cur = nxt
                    if len(loop) > len(bnd) + 2: ok = False; break
                if ok and len(loop) >= 3: loops.append(loop)
        if not loops: break
        newT = []
        for loop in loops:
            pts = V[loop]; c = pts.mean(0)
            _, _, vt = np.linalg.svd(pts - c); nrm = vt[2]
            P2 = _project(V[loop], nrm, c)
            tl = triangulate_with_holes([list(range(len(loop)))], P2)
            if not tl: continue
            newT.append(np.array(loop, int)[np.array(tl, int)])
        if not newT: break
        F = np.vstack([F] + newT)
        V, F = robust_weld(V, F, tol=5e-6 * diag)
    return V, F



# =============================== 그리드 절단 ===============================
def _plane_positions(lo, hi, step, eps):
    arr = np.arange(lo - eps, hi + eps, step)
    if arr[-1] < hi + eps:
        arr = np.append(arr, hi + eps)
    return arr


def grid_cut(V, F, sx, sy, sz, log=None, progress=None, fix_normals=True):
    """모델을 sx × sy × sz (mm) 격자로 자른다.
    반환: (pieces, axes)
      pieces : [{'ix','iy','iz','V','F'}, ...]
      axes   : {'x': xs, 'y': ys, 'z': zs}  (절단면 좌표 배열)
    """
    log = log or (lambda *a: None)
    mn = V.min(0); mx = V.max(0)
    eps = (mx - mn) * 1e-6 + 1e-9
    diag = np.linalg.norm(mx - mn) or 1.0
    tol = 1e-7 * diag

    xs = _plane_positions(mn[0], mx[0], sx, eps[0])
    ys = _plane_positions(mn[1], mx[1], sy, eps[1])
    zs = _plane_positions(mn[2], mx[2], sz, eps[2])
    nx, ny, nz = len(xs) - 1, len(ys) - 1, len(zs) - 1
    log(f"  격자 {nx} × {ny} × {nz} (X열 × Y행 × Z층) = 최대 {nx*ny*nz}칸")

    pieces = []
    done = 0
    total = max(1, nx * ny)
    for ix in range(nx):
        Vx, Fx = slab(V, F, 0, xs[ix], xs[ix + 1], tol, diag)
        if len(Fx) == 0:
            done += ny
            if progress: progress(done / total)
            continue
        for iy in range(ny):
            Vy, Fy = slab(Vx, Fx, 1, ys[iy], ys[iy + 1], tol, diag)
            done += 1
            if progress: progress(done / total)
            if len(Fy) == 0:
                continue
            for iz in range(nz):
                Vc, Fc = slab(Vy, Fy, 2, zs[iz], zs[iz + 1], tol, diag)
                if len(Fc) >= 4:
                    if fix_normals:
                        Fc = orient_faces(Vc, Fc)
                    pieces.append({"ix": ix, "iy": iy, "iz": iz, "V": Vc, "F": Fc})
    return pieces, {"x": xs, "y": ys, "z": zs}


# =============================== 법선 정리 ===============================
def signed_volume(V, F):
    """부호 있는 부피 (발산정리). 방향이 일관된 닫힌 메시에서만 의미가 있다."""
    if len(F) == 0:
        return 0.0
    T = V[F]
    return float(np.einsum("ij,ij->i", T[:, 0], np.cross(T[:, 1], T[:, 2])).sum() / 6.0)


def orient_faces(V, F):
    """절단/캡 과정에서 뒤집힌 삼각형의 방향을 이웃과 맞추고,
    덩어리마다 법선이 바깥을 향하도록 정리한다. (슬라이서 오류 방지)"""
    nf = len(F)
    if nf == 0:
        return F
    he = np.vstack([F[:, [0, 1]], F[:, [1, 2]], F[:, [2, 0]]])   # i = e*nf + f
    key = np.sort(he, axis=1)
    _, inv, _ = _unique_rows(key)
    ne = int(inv.max()) + 1
    order = np.argsort(inv, kind="stable")
    counts = np.bincount(inv, minlength=ne)
    starts = np.concatenate([[0], np.cumsum(counts)[:-1]])
    fwd = he[:, 0] < he[:, 1]

    visited = np.zeros(nf, bool)
    flip = np.zeros(nf, bool)
    from collections import deque
    for seed in range(nf):
        if visited[seed]:
            continue
        visited[seed] = True
        comp = [seed]
        dq = deque([seed])
        while dq:
            f = dq.popleft()
            for e in range(3):
                i = e * nf + f
                g_id = inv[i]
                if counts[g_id] != 2:            # 경계/비매니폴드 모서리는 건너뜀
                    continue
                s = starts[g_id]
                a, b = order[s], order[s + 1]
                j = b if a == i else a
                f2 = j % nf
                if visited[f2]:
                    continue
                visited[f2] = True
                flip[f2] = not (bool(fwd[j]) ^ bool(fwd[i]) ^ bool(flip[f]))
                comp.append(f2)
                dq.append(f2)
        # 덩어리 전체 법선 방향: 부피가 음수면 뒤집는다
        ci = np.array(comp, np.int64)
        Fc = F[ci].copy()
        fc = flip[ci]
        Fc[fc] = Fc[fc][:, ::-1]
        if signed_volume(V, Fc) < 0:
            flip[ci] = ~flip[ci]
    Fo = F.copy()
    Fo[flip] = Fo[flip][:, ::-1]
    return Fo


# =============================== 번호 매기기 ===============================
def number_pieces(pieces, axes):
    """아래층 -> 도면 위쪽 행(Y 큰 쪽) -> 왼쪽 열(X 작은 쪽) 순서로 1번부터 통번호를 매긴다."""
    order = sorted(range(len(pieces)),
                   key=lambda k: (pieces[k]["iz"], -pieces[k]["iy"], pieces[k]["ix"]))
    ny = len(axes["y"]) - 1
    layer_count = defaultdict(int)
    out = []
    for n, k in enumerate(order, start=1):
        p = dict(pieces[k])
        p["num"] = n
        p["layer"] = p["iz"] + 1                     # 1부터 시작하는 층 번호
        p["row"] = ny - p["iy"]                      # 도면 위쪽이 1행
        p["col"] = p["ix"] + 1                       # 왼쪽이 1열
        layer_count[p["iz"]] += 1
        p["in_layer"] = layer_count[p["iz"]]
        V = p["V"]
        p["bbox_min"] = V.min(0)
        p["bbox_max"] = V.max(0)
        out.append(p)
    return out


# =============================== 저장 ===============================
def save_pieces(pieces, out_dir, fmt="stl", log=None):
    """1.stl, 2.stl ... 처럼 숫자 파일명으로 저장한다."""
    log = log or (lambda *a: None)
    os.makedirs(out_dir, exist_ok=True)
    saved = []
    for p in pieces:
        name = f"{p['num']}.{fmt}"
        path = os.path.join(out_dir, name)
        save_mesh(path, p["V"], p["F"])
        rep = watertight_report(p["V"], p["F"])
        p["faces"] = int(len(p["F"]))
        p["watertight"] = bool(rep["watertight"])
        p["boundary_edges"] = int(rep["boundary_edges"])
        p["file"] = name
        saved.append(path)
        flag = "OK" if rep["watertight"] else f"경계{rep['boundary_edges']}"
        log(f"    {name:<10s} 층{p['layer']} 행{p['row']} 열{p['col']}  "
            f"면 {len(p['F']):>7,d}  [{flag}]")
    return saved


def write_index_csv(pieces, out_dir, filename="index.csv"):
    path = os.path.join(out_dir, filename)
    with open(path, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.writer(f)
        w.writerow(["number", "file", "layer", "row", "col",
                    "x_min_mm", "x_max_mm", "y_min_mm", "y_max_mm",
                    "z_min_mm", "z_max_mm", "faces", "watertight"])
        for p in pieces:
            a, b = p["bbox_min"], p["bbox_max"]
            w.writerow([p["num"], p.get("file", ""), p["layer"], p["row"], p["col"],
                        f"{a[0]:.3f}", f"{b[0]:.3f}", f"{a[1]:.3f}", f"{b[1]:.3f}",
                        f"{a[2]:.3f}", f"{b[2]:.3f}",
                        p.get("faces", len(p["F"])), int(p.get("watertight", False))])
    return path


# =============================== 위치 지도 ===============================
_PALETTE = ["#4E79A7", "#F28E2B", "#59A14F", "#E15759", "#B07AA1",
            "#76B7B2", "#EDC948", "#FF9DA7", "#9C755F", "#86BCB6"]


def _piece_color(num):
    return _PALETTE[(num - 1) % len(_PALETTE)]


def _hex_to_rgb(h):
    h = h.lstrip("#")
    return np.array([int(h[i:i + 2], 16) / 255.0 for i in (0, 2, 4)])


def _rasterize_layer(layer_pieces, x0, x1, y0, y1, px, density=12.0, max_samples=4_000_000):
    """층에 속한 조각들의 XY 단면(발자국)을 라벨 래스터로 만든다.
    반환: (H,W) int32, 값은 조각 번호(-1=빈칸)"""
    W = max(8, int(round((x1 - x0) / px)))
    H = max(8, int(round((y1 - y0) / px)))
    lab = np.full((H, W), -1, np.int32)
    px_area = ((x1 - x0) / W) * ((y1 - y0) / H)
    for p in layer_pieces:
        T = p["V"][p["F"]][:, :, :2]                       # (n,3,2) XY 투영
        a = T[:, 1] - T[:, 0]
        b = T[:, 2] - T[:, 0]
        area = np.abs(a[:, 0] * b[:, 1] - a[:, 1] * b[:, 0]) * 0.5
        ns = np.ceil(area / px_area * density).astype(np.int64)
        np.clip(ns, 0, 20000, out=ns)
        tot = int(ns.sum())
        if tot > max_samples and tot > 0:                  # 너무 많으면 균등 축소
            ns = np.maximum((ns * (max_samples / tot)).astype(np.int64), (area > 0).astype(np.int64))
            tot = int(ns.sum())
        pts = [T.reshape(-1, 2)]                           # 꼭짓점은 항상 포함
        if tot > 0:
            ti = np.repeat(np.arange(len(T)), ns)
            r1 = np.sqrt(np.random.random(tot))
            r2 = np.random.random(tot)
            w0 = (1.0 - r1)[:, None]
            w1 = (r1 * (1.0 - r2))[:, None]
            w2 = (r1 * r2)[:, None]
            pts.append(T[ti, 0] * w0 + T[ti, 1] * w1 + T[ti, 2] * w2)
        S = np.vstack(pts)
        cx = ((S[:, 0] - x0) / (x1 - x0) * W).astype(np.int64)
        cy = ((S[:, 1] - y0) / (y1 - y0) * H).astype(np.int64)
        np.clip(cx, 0, W - 1, out=cx)
        np.clip(cy, 0, H - 1, out=cy)
        m = np.zeros((H, W), bool)
        m[cy, cx] = True
        # 3x3 닫기(팽창->침식)로 샘플링 구멍 제거
        d = m.copy()
        for sy in (-1, 0, 1):
            for sx_ in (-1, 0, 1):
                d |= np.roll(np.roll(m, sy, 0), sx_, 1)
        e = d.copy()
        for sy in (-1, 0, 1):
            for sx_ in (-1, 0, 1):
                e &= np.roll(np.roll(d, sy, 0), sx_, 1)
        lab[e] = p["num"]
    return lab


def _label_xy(lab, num, x0, x1, y0, y1):
    """조각 번호 텍스트를 놓을 위치(도형 내부의 무게중심에 가까운 점)."""
    H, W = lab.shape
    ys, xs = np.nonzero(lab == num)
    if len(xs) == 0:
        return None
    cx, cy = xs.mean(), ys.mean()
    d = (xs - cx) ** 2 + (ys - cy) ** 2
    k = int(np.argmin(d))
    gx = x0 + (xs[k] + 0.5) / W * (x1 - x0)
    gy = y0 + (ys[k] + 0.5) / H * (y1 - y0)
    return gx, gy


def _lab_to_rgb(lab):
    H, W = lab.shape
    img = np.ones((H, W, 4))
    img[..., 3] = 0.0
    for num in np.unique(lab):
        if num < 0:
            continue
        m = lab == num
        img[m, :3] = _hex_to_rgb(_piece_color(int(num)))
        img[m, 3] = 0.85
    return img


def _draw_layer_axes(ax, lab, pieces_in_layer, axes, x0, x1, y0, y1, plt, label_size=13):
    ax.imshow(_lab_to_rgb(lab), origin="lower", extent=[x0, x1, y0, y1],
              interpolation="nearest", zorder=2)
    # 조각 외곽선
    for p in pieces_in_layer:
        m = (lab == p["num"]).astype(float)
        if m.max() == 0:
            continue
        ax.contour(np.linspace(x0, x1, lab.shape[1]),
                   np.linspace(y0, y1, lab.shape[0]),
                   m, levels=[0.5], colors="#1f2933", linewidths=1.1, zorder=3)
    # 절단 격자선
    for xv in axes["x"]:
        ax.axvline(xv, color="#9aa5b1", lw=0.7, ls=(0, (4, 3)), zorder=1)
    for yv in axes["y"]:
        ax.axhline(yv, color="#9aa5b1", lw=0.7, ls=(0, (4, 3)), zorder=1)
    # 번호
    for p in pieces_in_layer:
        pos = _label_xy(lab, p["num"], x0, x1, y0, y1)
        if pos is None:
            continue
        ax.text(pos[0], pos[1], str(p["num"]), ha="center", va="center",
                fontsize=label_size, fontweight="bold", color="#10171e", zorder=4,
                bbox=dict(boxstyle="round,pad=0.22", fc="white", ec="#1f2933",
                          lw=0.8, alpha=0.9))
    ax.set_xlim(x0, x1); ax.set_ylim(y0, y1)
    ax.set_aspect("equal")


def make_layer_maps(pieces, axes, out_dir, px_target=1100, dpi=150,
                    log=None, make_overview=True):
    """각 층을 위에서 내려다본 평면도 PNG + 전체 요약 PNG 를 만든다."""
    log = log or (lambda *a: None)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    allmin = np.min([p["bbox_min"] for p in pieces], axis=0)
    allmax = np.max([p["bbox_max"] for p in pieces], axis=0)
    pad = 0.04 * max(allmax[0] - allmin[0], allmax[1] - allmin[1], 1e-6)
    x0, x1 = allmin[0] - pad, allmax[0] + pad
    y0, y1 = allmin[1] - pad, allmax[1] + pad
    px = max(x1 - x0, y1 - y0) / px_target

    zs = axes["z"]
    layers = sorted({p["iz"] for p in pieces})
    rasters = {}
    paths = []
    for iz in layers:
        lp = [p for p in pieces if p["iz"] == iz]
        lab = _rasterize_layer(lp, x0, x1, y0, y1, px)
        rasters[iz] = lab
        w_in = 10.0
        h_in = w_in * (y1 - y0) / (x1 - x0)
        h_in = float(np.clip(h_in, 3.0, 16.0))
        fig, ax = plt.subplots(figsize=(w_in, h_in + 0.9))
        _draw_layer_axes(ax, lab, lp, axes, x0, x1, y0, y1, plt)
        nums = [p["num"] for p in lp]
        ax.set_title(
            f"LAYER {iz + 1} of {len(zs) - 1}   |   Z = {max(zs[iz], allmin[2]):.1f} – "
            f"{min(zs[iz + 1], allmax[2]):.1f} mm   |   {len(lp)} pieces "
            f"(#{min(nums)}–#{max(nums)})",
            fontsize=13, fontweight="bold", pad=12)
        ax.set_xlabel("X (mm)  →", fontsize=10)
        ax.set_ylabel("Y (mm)  ↑", fontsize=10)
        ax.tick_params(labelsize=8)
        for s in ax.spines.values():
            s.set_color("#cbd2d9")
        fig.text(0.01, 0.01, "Top view (looking down the −Z direction). "
                             "Dashed lines = cut planes.",
                 fontsize=8, color="#616e7c")
        fig.tight_layout()
        path = os.path.join(out_dir, f"layer_{iz + 1:02d}.png")
        fig.savefig(path, dpi=dpi, facecolor="white")
        plt.close(fig)
        paths.append(path)
        log(f"    지도 저장: {os.path.basename(path)}  (조각 {len(lp)}개)")

    if make_overview and len(layers) > 1:
        ncol = min(3, len(layers))
        nrow = int(math.ceil(len(layers) / ncol))
        cw = 5.0
        ch = cw * (y1 - y0) / (x1 - x0)
        ch = float(np.clip(ch, 2.2, 8.0))
        fig, axs = plt.subplots(nrow, ncol, figsize=(cw * ncol, (ch + 0.7) * nrow),
                                squeeze=False)
        for k, iz in enumerate(layers):
            ax = axs[k // ncol][k % ncol]
            lp = [p for p in pieces if p["iz"] == iz]
            _draw_layer_axes(ax, rasters[iz], lp, axes, x0, x1, y0, y1, plt, label_size=9)
            ax.set_title(f"Layer {iz + 1}  (Z {max(zs[iz], allmin[2]):.0f}–"
                         f"{min(zs[iz + 1], allmax[2]):.0f} mm)",
                         fontsize=10, fontweight="bold")
            ax.tick_params(labelsize=7)
        for k in range(len(layers), nrow * ncol):
            axs[k // ncol][k % ncol].axis("off")
        fig.suptitle("ASSEMBLY MAP — all layers, top view (bottom layer first)",
                     fontsize=14, fontweight="bold")
        fig.tight_layout(rect=(0, 0, 1, 0.97))
        path = os.path.join(out_dir, "layers_overview.png")
        fig.savefig(path, dpi=dpi, facecolor="white")
        plt.close(fig)
        paths.append(path)
        log(f"    지도 저장: {os.path.basename(path)}  (전체 {len(layers)}개 층)")
    return paths


def make_exploded_map(pieces, axes, out_dir, explode=0.25, dpi=150,
                      max_tris_total=150_000, log=None):
    """조각을 벌려 놓은 3D 분해 조립도 PNG.
    층 방향(Z)으로 크게 벌려서 각 층이 확실히 구분되도록 그린다."""
    log = log or (lambda *a: None)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection

    allmin = np.min([p["bbox_min"] for p in pieces], axis=0)
    allmax = np.max([p["bbox_max"] for p in pieces], axis=0)
    ctr = (allmin + allmax) / 2.0
    size = np.maximum(allmax - allmin, 1e-9)
    step = np.array([axes["x"][1] - axes["x"][0],
                     axes["y"][1] - axes["y"][0],
                     axes["z"][1] - axes["z"][0]])

    # 조각별 이동량: 평면 방향은 조금, 층 방향은 크게 벌린다.
    gap = np.array([explode * 0.55, explode * 0.55, max(explode, 0.35) * 1.6])
    shifts = {}
    for p in pieces:
        cell_ctr = (p["bbox_min"] + p["bbox_max"]) / 2.0
        direction = (cell_ctr - ctr) / step
        s = direction * step * gap
        s[2] = (p["iz"] - (len(axes["z"]) - 2) / 2.0) * step[2] * gap[2]
        shifts[p["num"]] = s

    # 삼각형 총량 예산 배분 (희박한 줄무늬가 생기지 않도록 넉넉하게)
    tot_f = sum(len(p["F"]) for p in pieces)
    budget = max(1.0, min(1.0, max_tris_total / max(tot_f, 1)))

    fig = plt.figure(figsize=(11, 10))
    try:
        ax = fig.add_subplot(111, projection="3d", computed_zorder=False)
    except Exception:
        ax = fig.add_subplot(111, projection="3d")
    elev, azim = 24.0, -60.0
    ve = math.radians(elev); va = math.radians(azim)
    view = np.array([math.cos(ve) * math.cos(va), math.cos(ve) * math.sin(va),
                     math.sin(ve)])

    def depth(p):
        c = (p["bbox_min"] + p["bbox_max"]) / 2.0 + shifts[p["num"]]
        return float(np.dot(c, view))

    ordered = sorted(pieces, key=depth)          # 뒤쪽부터 그린다
    labels = []
    for k, p in enumerate(ordered):
        Vs = p["V"] + shifts[p["num"]]
        Fp = p["F"]
        st = max(1, int(round(1.0 / budget)))
        polys = Vs[Fp[::st]] if st > 1 else Vs[Fp]
        col = Poly3DCollection(polys, facecolors=_piece_color(p["num"]),
                               alpha=1.0, shade=True)
        col.set_linewidth(0)
        col.set_zsort("average")
        col.set_zorder(k + 2)
        ax.add_collection3d(col)
        c = (Vs.min(0) + Vs.max(0)) / 2.0
        labels.append((k, c, p["num"]))
    for k, c, num in labels:
        ax.text(c[0], c[1], c[2], str(num), fontsize=10, fontweight="bold",
                ha="center", va="center", color="#10171e", zorder=len(ordered) + 10,
                bbox=dict(boxstyle="round,pad=0.18", fc="white", ec="#1f2933",
                          lw=0.7, alpha=0.92))

    lo = np.min([p["bbox_min"] + shifts[p["num"]] for p in pieces], axis=0)
    hi = np.max([p["bbox_max"] + shifts[p["num"]] for p in pieces], axis=0)
    pad = (hi - lo) * 0.06 + 1e-6
    lo -= pad; hi += pad
    ax.set_xlim(lo[0], hi[0]); ax.set_ylim(lo[1], hi[1]); ax.set_zlim(lo[2], hi[2])
    try:
        ax.set_box_aspect(tuple(hi - lo))
    except Exception:
        pass
    ax.set_xlabel("X (mm)", fontsize=9); ax.set_ylabel("Y (mm)", fontsize=9)
    ax.set_zlabel("Z (mm)", fontsize=9)
    ax.tick_params(labelsize=7)
    for pane in (ax.xaxis, ax.yaxis, ax.zaxis):
        pane.pane.set_facecolor("white")
        pane.pane.set_edgecolor("#e4e7eb")
        pane._axinfo["grid"]["color"] = "#eef1f4"
    n_layers = len({p["iz"] for p in pieces})
    ax.set_title(f"EXPLODED ASSEMBLY VIEW — {len(pieces)} pieces in {n_layers} layer(s)",
                 fontsize=13, fontweight="bold")
    fig.text(0.5, 0.02, "Layers are pulled apart vertically. "
                        "Piece numbers run bottom layer first.",
             ha="center", fontsize=8, color="#616e7c")
    ax.view_init(elev=elev, azim=azim)
    fig.tight_layout()
    path = os.path.join(out_dir, "exploded_3d.png")
    fig.savefig(path, dpi=dpi, facecolor="white")
    plt.close(fig)
    log(f"    지도 저장: {os.path.basename(path)}")
    return path


# =============================== 전체 실행 ===============================
def run_split(input_path, output_dir,
              cut_x=200.0, cut_y=200.0, cut_z=200.0,
              output_format="stl",
              repair_source_holes=True, fix_normals=True,
              layer_maps=True, exploded_map=True, index_csv=True,
              make_zip=False, explode=0.25,
              log=print, progress=None):
    """모든 단계를 순서대로 실행한다. GUI/CLI 공통 진입점."""
    log = log or (lambda *a: None)
    progress = progress or (lambda f: None)

    if not os.path.exists(input_path):
        raise FileNotFoundError(f"입력 파일을 찾을 수 없습니다: {input_path}")
    os.makedirs(output_dir, exist_ok=True)

    progress(0.02)
    log(f"[1/5] 불러오는 중: {os.path.basename(input_path)}")
    V, F = load_mesh(input_path)
    rep = watertight_report(V, F)
    log(f"      정점 {len(V):,} / 면 {len(F):,} | 원본 경계모서리 {rep['boundary_edges']}")

    if repair_source_holes and rep["boundary_edges"] > 0:
        log("[2/5] 원본 구멍 메우는 중...")
        V, F = fill_holes(V, F)
        r2 = watertight_report(V, F)
        log(f"      → 경계모서리 {r2['boundary_edges']}, 비매니폴드 {r2['nonmanifold']}")
    else:
        log("[2/5] 구멍 메우기 건너뜀")
    progress(0.10)

    log(f"[3/5] 절단 중: {cut_x} × {cut_y} × {cut_z} mm")
    raw, axes = grid_cut(V, F, cut_x, cut_y, cut_z, log=log, fix_normals=fix_normals,
                         progress=lambda f: progress(0.10 + 0.55 * f))
    if not raw:
        raise RuntimeError("조각이 하나도 생성되지 않았습니다. 절단 크기를 확인하세요.")
    pieces = number_pieces(raw, axes)
    n_layers = len({p["iz"] for p in pieces})
    log(f"      조각 {len(pieces)}개 / 층 {n_layers}개")
    progress(0.66)

    log(f"[4/5] 저장 중 → {output_dir}")
    saved = save_pieces(pieces, output_dir, fmt=output_format, log=log)
    progress(0.80)

    maps = []
    csv_path = None
    if index_csv:
        csv_path = write_index_csv(pieces, output_dir)
        log(f"      목록 저장: {os.path.basename(csv_path)}")
    log("[5/5] 위치 지도 만드는 중...")
    if layer_maps:
        try:
            maps += make_layer_maps(pieces, axes, output_dir, log=log)
        except Exception as ex:
            log(f"      (평면도 생략: {ex})")
    progress(0.92)
    if exploded_map:
        try:
            maps.append(make_exploded_map(pieces, axes, output_dir,
                                          explode=explode, log=log))
        except Exception as ex:
            log(f"      (3D 조립도 생략: {ex})")

    zip_path = None
    if make_zip:
        zip_path = os.path.join(output_dir, "pieces.zip")
        with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as z:
            for p in saved:
                z.write(p, os.path.basename(p))
            for m in maps:
                z.write(m, os.path.basename(m))
            if csv_path:
                z.write(csv_path, os.path.basename(csv_path))
        log(f"      ZIP 생성: {os.path.basename(zip_path)}")

    progress(1.0)
    ok = sum(1 for p in pieces if p.get("watertight"))
    log(f"완료! 조각 {len(pieces)}개 (방수 {ok}개), 층 {n_layers}개, 지도 {len(maps)}장")
    return {"pieces": pieces, "axes": axes, "files": saved,
            "maps": maps, "csv": csv_path, "zip": zip_path,
            "n_layers": n_layers}


# =============================== CLI ===============================
def main(argv=None):
    import argparse
    ap = argparse.ArgumentParser(
        description="3D 모델을 격자로 잘라 숫자 파일명으로 저장하고 위치 지도를 만듭니다.")
    ap.add_argument("input", help="입력 파일 (.stl / .obj)")
    ap.add_argument("-o", "--output", default="pieces", help="출력 폴더")
    ap.add_argument("-x", type=float, default=200.0, help="X 절단 크기 (mm)")
    ap.add_argument("-y", type=float, default=200.0, help="Y 절단 크기 (mm)")
    ap.add_argument("-z", type=float, default=200.0, help="Z 절단 크기 (mm)")
    ap.add_argument("-f", "--format", default="stl", choices=["stl", "obj"])
    ap.add_argument("--no-repair", action="store_true", help="원본 구멍 메우기 끄기")
    ap.add_argument("--no-fix-normals", action="store_true",
                    help="조각 법선 방향 자동 정리 끄기")
    ap.add_argument("--no-maps", action="store_true", help="지도 생성 끄기")
    ap.add_argument("--no-exploded", action="store_true", help="3D 조립도 끄기")
    ap.add_argument("--explode", type=float, default=0.25, help="3D 조립도 벌림 정도")
    ap.add_argument("--zip", action="store_true", help="결과를 ZIP 으로 묶기")
    a = ap.parse_args(argv)
    run_split(a.input, a.output, a.x, a.y, a.z,
              output_format=a.format,
              repair_source_holes=not a.no_repair,
              fix_normals=not a.no_fix_normals,
              layer_maps=not a.no_maps,
              exploded_map=not (a.no_maps or a.no_exploded),
              make_zip=a.zip, explode=a.explode)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
