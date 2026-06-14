# 필요한 외부 패키지 설치 및 임포트
!pip install numpy scipy matplotlib

import os
import sys
import struct
from collections import deque, defaultdict
import numpy as np
from scipy.spatial import cKDTree
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
# @title 3D 모델 지정 크기(mm) 절단 실행기
# @markdown 아래 항목에 모델 경로와 자르고자 하는 각 조각의 최대 크기(mm)를 입력하세요. 
# @markdown (예: 프린터 베드가 200mm라면 여유를 두어 180~190으로 설정)

# ----------------------------------------------------------------------------
# 1. 사용자 설정 (Colab Form)
# ----------------------------------------------------------------------------
input_file_path = "/content/drive/MyDrive/3D파일(박물관)/MOR555_SCALE.stl" # @param {type:"string"}
output_directory_path = "/content/drive/MyDrive/3D파일(박물관)/티라노" # @param {type:"string"}
cut_size_x_mm = 400 # @param {type:"number"}
cut_size_y_mm = 400 # @param {type:"number"}
cut_size_z_mm = 400 # @param {type:"number"}
download_result_zip = True # @param {type:"boolean"}
show_preview = True # @param {type:"boolean"}


# ----------------------------------------------------------------------------
# 2. 핵심 함수 정의 (단위 변환 및 절단 알고리즘)
# ----------------------------------------------------------------------------

def load_mesh(path):
    ext = os.path.splitext(path)[1].lower()
    if ext == ".stl": return _load_stl(path)
    if ext == ".obj": return _load_obj(path)
    raise ValueError(f"지원하지 않는 포맷: {ext} (stl/obj 지원)")

def save_mesh(path, V, F):
    ext = os.path.splitext(path)[1].lower()
    if ext == ".stl": return _save_stl_binary(path, V, F)
    if ext == ".obj": return _save_obj(path, V, F)
    raise ValueError(f"지원하지 않는 포맷: {ext}")

def _load_stl(path):
    with open(path, "rb") as f:
        head = f.read(5)
        f.seek(0)
        data = f.read()
    if head == b"solid" and b"facet" in data[:2048]:
        return _load_stl_ascii(data.decode("utf-8", "ignore"))
    return _load_stl_binary(data)

def _load_stl_binary(data):
    n = struct.unpack("<I", data[80:84])[0]
    tris = np.zeros((n, 3, 3), dtype=np.float64)
    off = 84
    for i in range(n):
        vals = struct.unpack("<12fH", data[off:off + 50])
        tris[i] = np.array(vals[3:12]).reshape(3, 3)
        off += 50
    V = tris.reshape(-1, 3)
    F = np.arange(len(V)).reshape(-1, 3)
    return weld(V, F)

def _load_stl_ascii(text):
    verts = []
    for line in text.splitlines():
        line = line.strip()
        if line.startswith("vertex"):
            verts.append([float(x) for x in line.split()[1:4]])
    V = np.array(verts, dtype=np.float64)
    F = np.arange(len(V)).reshape(-1, 3)
    return weld(V, F)

def _load_obj(path):
    verts, faces = [], []
    with open(path) as f:
        for line in f:
            if line.startswith("v "):
                verts.append([float(x) for x in line.split()[1:4]])
            elif line.startswith("f "):
                idx = [int(p.split("/")[0]) - 1 for p in line.split()[1:]]
                for k in range(1, len(idx) - 1):
                    faces.append([idx[0], idx[k], idx[k + 1]])
    return np.array(verts, dtype=np.float64), np.array(faces, dtype=np.int64)

def _save_stl_binary(path, V, F):
    tris = V[F]
    e1 = tris[:, 1] - tris[:, 0]
    e2 = tris[:, 2] - tris[:, 0]
    nor = np.cross(e1, e2)
    ln = np.linalg.norm(nor, axis=1, keepdims=True)
    ln[ln == 0] = 1.0
    nor = nor / ln
    with open(path, "wb") as f:
        f.write(b"\0" * 80)
        f.write(struct.pack("<I", len(F)))
        for i in range(len(F)):
            f.write(struct.pack("<3f", *nor[i]))
            for v in tris[i]:
                f.write(struct.pack("<3f", *v))
            f.write(struct.pack("<H", 0))

def _save_obj(path, V, F):
    with open(path, "w") as f:
        for v in V:
            f.write(f"v {v[0]} {v[1]} {v[2]}\n")
        for tri in F:
            f.write(f"f {tri[0] + 1} {tri[1] + 1} {tri[2] + 1}\n")

def weld(V, F, tol=1e-7):
    if len(V) == 0: return V, F
    diag = np.linalg.norm(V.max(0) - V.min(0)) or 1.0
    q = np.round(V / (tol * diag)).astype(np.int64)
    _, inv, idx = _unique_rows(q)
    Vn = V[idx]
    Fn = inv[F]
    good = (Fn[:, 0] != Fn[:, 1]) & (Fn[:, 1] != Fn[:, 2]) & (Fn[:, 0] != Fn[:, 2])
    return Vn, Fn[good]

def _unique_rows(a):
    order = np.lexsort(a.T[::-1])
    a_sorted = a[order]
    diff = np.ones(len(a), bool)
    diff[1:] = np.any(a_sorted[1:] != a_sorted[:-1], axis=1)
    uid = np.cumsum(diff) - 1
    inv = np.empty(len(a), np.int64)
    inv[order] = uid
    first_idx = order[diff]
    return None, inv, first_idx

def watertight_report(V, F):
    if len(F) == 0:
        return {"watertight": False, "faces": 0, "boundary_edges": 0, "nonmanifold": 0}
    e = np.sort(np.vstack([F[:, [0, 1]], F[:, [1, 2]], F[:, [2, 0]]]), axis=1)
    _, inv, _ = _unique_rows(e)
    counts = np.bincount(inv)
    boundary = int(np.sum(counts == 1))
    nonman = int(np.sum(counts > 2))
    return {
        "watertight": boundary == 0 and nonman == 0,
        "faces": len(F),
        "boundary_edges": boundary,
        "nonmanifold": nonman,
    }

def reorient(V, F):
    if len(F) == 0: return V, F
    F = F.copy()
    nf = len(F)
    he = np.vstack([F[:, [0, 1]], F[:, [1, 2]], F[:, [2, 0]]])
    face_of = np.tile(np.arange(nf), 3)
    base_sign = (he[:, 0] < he[:, 1])
    key = np.sort(he, axis=1)
    _, inv, _ = _unique_rows(key)
    order = np.argsort(inv, kind="stable")
    inv_s = inv[order]
    bounds = np.flatnonzero(np.diff(inv_s)) + 1
    groups = np.split(order, bounds)

    adj = [[] for _ in range(nf)]
    for g in groups:
        if len(g) == 2:
            i, j = g
            fi, fj = face_of[i], face_of[j]
            adj[fi].append((fj, i, j))
            adj[fj].append((fi, j, i))

    flip = np.zeros(nf, bool)
    visited = np.zeros(nf, bool)
    for seed in range(nf):
        if visited[seed]: continue
        visited[seed] = True
        dq = deque([seed])
        while dq:
            f = dq.popleft()
            for (g, hi, hj) in adj[f]:
                if visited[g]: continue
                visited[g] = True
                eff_f = base_sign[hi] ^ flip[f]
                eff_g = base_sign[hj] ^ flip[g]
                if eff_f == eff_g:
                    flip[g] = True
                dq.append(g)
    F[flip] = F[flip][:, ::-1]

    rep = watertight_report(V, F)
    if rep["boundary_edges"] == 0:
        vol = np.einsum("ij,ij->i", V[F[:, 0]],
                        np.cross(V[F[:, 1]], V[F[:, 2]])).sum() / 6.0
        if vol < 0:
            F = F[:, ::-1]
    return V, F

def clip_cap(V, F, p0, n, tol=None):
    if len(F) == 0: return V, F
    n = np.asarray(n, float)
    n = n / (np.linalg.norm(n) or 1.0)
    p0 = np.asarray(p0, float)
    if tol is None:
        diag = np.linalg.norm(V.max(0) - V.min(0)) or 1.0
        tol = 1e-9 * diag

    d = (V - p0) @ n
    d[np.abs(d) < tol] = 0.0

    out_tris = []
    for tri in F:
        verts = V[tri]
        dist = d[tri]
        poly = _clip_triangle(verts, dist)
        for k in range(1, len(poly) - 1):
            out_tris.append([poly[0], poly[k], poly[k + 1]])

    if not out_tris:
        return np.zeros((0, 3)), np.zeros((0, 3), int)

    arr = np.array(out_tris)
    Vc = arr.reshape(-1, 3)
    Fc = np.arange(len(Vc)).reshape(-1, 3)
    Vc, Fc = weld(Vc, Fc, tol=1e-6)
    Vc, Fc = reorient(Vc, Fc)

    dc = (Vc - p0) @ n
    on_plane = np.abs(dc) < (tol * 1e3 + 1e-9)

    he = np.vstack([Fc[:, [0, 1]], Fc[:, [1, 2]], Fc[:, [2, 0]]])
    key = np.sort(he, axis=1)
    _, inv, _ = _unique_rows(key)
    counts = np.bincount(inv, minlength=inv.max() + 1)
    is_boundary = counts[inv] == 1
    sel = is_boundary & on_plane[he[:, 0]] & on_plane[he[:, 1]]
    cap_edges = he[sel]

    cap_tris = _cap_loops(Vc, cap_edges, n, p0)
    if len(cap_tris):
        Fc = np.vstack([Fc, cap_tris])

    return weld(Vc, Fc, tol=1e-6)

def _clip_triangle(verts, dist):
    out = []
    for i in range(3):
        cur, nxt = verts[i], verts[(i + 1) % 3]
        dc, dn = dist[i], dist[(i + 1) % 3]
        if dc >= 0: out.append(cur)
        if (dc < 0) != (dn < 0) and dc != dn:
            t = dc / (dc - dn)
            out.append(cur + t * (nxt - cur))
    return out

_CUR_P2 = None

def _project(V, n, p0):
    n = np.asarray(n, float)
    n = n / (np.linalg.norm(n) or 1.0)
    u = np.cross(n, [1., 0, 0])
    if np.linalg.norm(u) < 1e-6: u = np.cross(n, [0, 1., 0])
    u /= np.linalg.norm(u)
    w = np.cross(n, u)
    P2 = np.column_stack([(V - p0) @ u, (V - p0) @ w])
    return P2, u, w

def _signed_area(poly2):
    x = poly2[:, 0]; y = poly2[:, 1]
    return 0.5 * np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y)

def _assemble_loops(edges):
    P2 = _CUR_P2
    out = defaultdict(list)
    for k, (a, b) in enumerate(edges):
        out[a].append([b, 1])

    def ang(a, b):
        d = P2[b] - P2[a]
        return np.arctan2(d[1], d[0])

    loops = []
    remaining = len(edges)
    while remaining > 0:
        start = None
        for v, lst in out.items():
            if any(av for _, av in lst):
                start = v
                break
        if start is None: break
        loop = [start]
        cur = start
        prev = None
        for e in out[cur]:
            if e[1]:
                e[1] = 0; nxt = e[0]; break
        remaining -= 1
        loop.append(nxt)
        prev = cur; cur = nxt
        while cur != start:
            inc = ang(prev, cur)
            best = None; bestkey = None
            for e in out[cur]:
                if not e[1]: continue
                a = ang(cur, e[0])
                turn = (a - inc) % (2 * np.pi)
                if bestkey is None or turn < bestkey:
                    bestkey = turn; best = e
            if best is None: break
            best[1] = 0; remaining -= 1
            prev = cur; cur = best[0]
            loop.append(cur)
            if len(loop) > len(edges) + 2: break
        if loop[0] == loop[-1] and len(loop) > 3:
            loops.append(loop[:-1])
        elif len(loop) >= 3:
            loops.append(loop[:-1] if loop[0] == loop[-1] else loop)
    return loops

def _point_in_tri(p, a, b, c):
    d1 = (b[0] - a[0]) * (p[1] - a[1]) - (b[1] - a[1]) * (p[0] - a[0])
    d2 = (c[0] - b[0]) * (p[1] - b[1]) - (c[1] - b[1]) * (p[0] - b[0])
    d3 = (a[0] - c[0]) * (p[1] - c[1]) - (a[1] - c[1]) * (p[0] - c[0])
    neg = (d1 < 0) or (d2 < 0) or (d3 < 0)
    pos = (d1 > 0) or (d2 > 0) or (d3 > 0)
    return not (neg and pos)

def _ear_clip(poly_idx, P2):
    idx = list(poly_idx)
    if len(idx) < 3: return []
    if _signed_area(P2[idx]) < 0: idx = idx[::-1]
    pts = P2; tris = []

    def cross(a, b, c):
        return (pts[b, 0] - pts[a, 0]) * (pts[c, 1] - pts[a, 1]) - \
               (pts[b, 1] - pts[a, 1]) * (pts[c, 0] - pts[a, 0])

    guard = 0
    while len(idx) > 3 and guard < 5 * len(poly_idx) + 10000:
        guard += 1; m = len(idx); found = False
        for i in range(m):
            a, b, c = idx[(i - 1) % m], idx[i], idx[(i + 1) % m]
            if cross(a, b, c) <= 1e-14: continue
            pa, pb, pc = pts[a], pts[b], pts[c]; ok = True
            for j in idx:
                if j in (a, b, c): continue
                if _point_in_tri(pts[j], pa, pb, pc):
                    ok = False; break
            if ok:
                tris.append((a, b, c)); del idx[i]; found = True; break
        if not found:
            best = None; bestv = -1
            for i in range(m):
                a, b, c = idx[(i - 1) % m], idx[i], idx[(i + 1) % m]
                cv = cross(a, b, c)
                if cv > bestv: bestv = cv; best = i
            if best is None or bestv <= 0: break
            i = best
            a, b, c = idx[(i - 1) % m], idx[i], idx[(i + 1) % m]
            tris.append((a, b, c)); del idx[i]
    if len(idx) == 3: tris.append((idx[0], idx[1], idx[2]))
    return tris

def _find_visible(hole, outer, P2):
    hi = max(range(len(hole)), key=lambda i: P2[hole[i], 0])
    M = P2[hole[hi]]
    bestI = None; bestx = np.inf; edge = None; m = len(outer)
    for i in range(m):
        a = P2[outer[i]]; b = P2[outer[(i + 1) % m]]
        ya, yb = a[1], b[1]
        if (ya > M[1]) == (yb > M[1]): continue
        t = (M[1] - ya) / (yb - ya)
        xi = a[0] + t * (b[0] - a[0])
        if xi >= M[0] - 1e-12 and xi < bestx:
            bestx = xi; bestI = np.array([xi, M[1]])
            edge = (outer[i], outer[(i + 1) % m])
    if edge is None:
        return min(range(m), key=lambda i: np.sum((P2[outer[i]] - M) ** 2)), hi
    P = edge[0] if P2[edge[0], 0] > P2[edge[1], 0] else edge[1]
    Pp = P2[P]; cand = P; bestang = -np.inf
    for v in outer:
        if v == P: continue
        pv = P2[v]
        if _point_in_tri(pv, M, bestI, Pp):
            d = pv - M
            ang = d[0] / (np.linalg.norm(d) + 1e-30)
            if ang > bestang: bestang = ang; cand = v
    vis_idx = outer.index(cand)
    return vis_idx, hi

def _bridge(outer, hole, P2):
    vis_idx, hi = _find_visible(hole, outer, P2)
    hrot = hole[hi:] + hole[:hi]
    new = outer[:vis_idx + 1] + hrot + [hrot[0]] + outer[vis_idx:]
    return new

def _cap_loops(V, cap_edges, n, p0):
    if len(cap_edges) == 0: return np.zeros((0, 3), int)
    P2, u, w = _project(V, n, p0)
    global _CUR_P2
    _CUR_P2 = P2
    loops = _assemble_loops([tuple(e) for e in cap_edges.tolist()])
    if not loops: return np.zeros((0, 3), int)
    areas = [_signed_area(P2[l]) for l in loops]

    def pip(p, poly):
        x, y = p; inside = False; m = len(poly); j = m - 1
        for i in range(m):
            xi, yi = poly[i]; xj, yj = poly[j]
            if ((yi > y) != (yj > y)) and (x < (xj - xi) * (y - yi) / (yj - yi + 1e-30) + xi):
                inside = not inside
            j = i
        return inside

    reps = [P2[l[0]] for l in loops]
    depth = [0] * len(loops)
    for i in range(len(loops)):
        for j in range(len(loops)):
            if i == j: continue
            if pip(reps[i], P2[loops[j]]): depth[i] += 1
    outers = [i for i in range(len(loops)) if depth[i] % 2 == 0]
    holes = [i for i in range(len(loops)) if depth[i] % 2 == 1]
    L = [list(l) for l in loops]
    for i in outers:
        if areas[i] < 0: L[i] = L[i][::-1]
    for i in holes:
        if areas[i] > 0: L[i] = L[i][::-1]
    tris = []
    for oi in outers:
        opoly = P2[L[oi]]
        mine = [hi for hi in holes if pip(reps[hi], opoly) and (depth[hi] == depth[oi] + 1)]
        merged = list(L[oi])
        for hi in sorted(mine, key=lambda h: -max(P2[L[h]][:, 0])):
            merged = _bridge(merged, L[hi], P2)
        tris += _ear_clip(merged, P2)
    if not tris: return np.zeros((0, 3), int)
    T = np.array(tris, int)
    nrm = np.cross(V[T[:, 1]] - V[T[:, 0]], V[T[:, 2]] - V[T[:, 0]])
    flip = (nrm @ n) > 0
    T[flip] = T[flip][:, ::-1]
    return T

def robust_weld(V, F, tol):
    if len(V) == 0: return V, F
    tree = cKDTree(V)
    pairs = tree.query_pairs(r=tol, output_type='ndarray')
    parent = np.arange(len(V))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for i, j in pairs:
        ri, rj = find(i), find(j)
        if ri != rj:
            parent[max(ri, rj)] = min(ri, rj)
    root = np.array([find(i) for i in range(len(V))])
    uniq, inv = np.unique(root, return_inverse=True)
    Vn = V[uniq]; Fn = inv[F]
    good = (Fn[:, 0] != Fn[:, 1]) & (Fn[:, 1] != Fn[:, 2]) & (Fn[:, 0] != Fn[:, 2])
    return Vn, Fn[good]

def _dedup_degenerate(V, F, area_tol):
    if len(F) == 0: return F
    a = V[F[:, 0]]; b = V[F[:, 1]]; c = V[F[:, 2]]
    area = 0.5 * np.linalg.norm(np.cross(b - a, c - a), axis=1)
    F = F[area > area_tol]
    if len(F) == 0: return F
    st = np.sort(F, axis=1)
    _, inv, first = _unique_rows(st)
    cnt = np.bincount(inv)
    keep = first[cnt[np.arange(len(first))] % 2 == 1]
    return F[keep]

def heal_tjunctions(V, F, tol=None):
    if len(F) == 0: return V, F
    if tol is None: tol = 1e-7 * (np.linalg.norm(V.max(0) - V.min(0)) or 1.0)
    F = F.tolist()
    for _ in range(200):
        Fa = np.array(F)
        he = np.vstack([Fa[:, [0, 1]], Fa[:, [1, 2]], Fa[:, [2, 0]]])
        key = np.sort(he, axis=1)
        _, inv, _ = _unique_rows(key)
        cnt = np.bincount(inv)
        bmask = cnt[inv] == 1; bnd = he[bmask]
        if len(bnd) == 0: break
        bverts = np.unique(bnd.reshape(-1)); split = None
        for (a, b) in bnd:
            pa, pb = V[a], V[b]; ab = pb - pa; L2 = ab @ ab
            if L2 < tol * tol: continue
            for c in bverts:
                if c == a or c == b: continue
                pc = V[c]; t = ((pc - pa) @ ab) / L2
                if t <= 1e-6 or t >= 1 - 1e-6: continue
                proj = pa + t * ab
                if np.linalg.norm(pc - proj) < tol * 50:
                    split = (a, b, c); break
            if split: break
        if split is None: break
        a, b, c = split; newF = []; done = False
        for tri in F:
            if not done and a in tri and b in tri:
                third = [v for v in tri if v != a and v != b]
                if len(third) != 1:
                    newF.append(tri); continue
                seq = tri; res = []; m = len(seq)
                for i in range(m):
                    uu = seq[i]; ww = seq[(i + 1) % m]
                    res.append(uu)
                    if (uu == a and ww == b) or (uu == b and ww == a):
                        res.append(c)
                if len(res) == 4:
                    newF.append([res[0], res[1], res[2]])
                    newF.append([res[0], res[2], res[3]])
                    done = True
                else: newF.append(tri)
            else: newF.append(tri)
        if not done: break
        F = newF
    return V, np.array(F, int)

def finalize(V, F, tol=None):
    if len(F) == 0: return V, F
    diag = np.linalg.norm(V.max(0) - V.min(0)) or 1.0
    htol = 1e-6 * diag
    for _ in range(8):
        V, F = robust_weld(V, F, tol=5e-6 * diag)
        F = _dedup_degenerate(V, F, 1e-10 * diag * diag)
        V, F = heal_tjunctions(V, F, tol=htol)
        if watertight_report(V, F)['watertight']: break
    V, F = reorient(V, F)
    return V, F

def clip_slab(V, F, axis, lo, hi):
    nplus = np.zeros(3); nplus[axis] = 1.0; p_lo = np.zeros(3); p_lo[axis] = lo
    V, F = clip_cap(V, F, p_lo, nplus)
    if len(F) == 0: return V, F
    p_hi = np.zeros(3); p_hi[axis] = hi
    V, F = clip_cap(V, F, p_hi, -nplus)
    return V, F

def grid_cut_by_size(V, F, sx, sy, sz, do_finalize=True):
    """지정된 크기(mm) 간격으로 메쉬를 자릅니다."""
    mn = V.min(0); mx = V.max(0)
    eps = (mx - mn) * 1e-6 + 1e-9

    def make_planes(start, end, step, e):
        arr = np.arange(start - e, end + e, step)
        if arr[-1] < end + e:
            arr = np.append(arr, end + e)
        return arr

    xs = make_planes(mn[0], mx[0], sx, eps[0])
    ys = make_planes(mn[1], mx[1], sy, eps[1])
    zs = make_planes(mn[2], mx[2], sz, eps[2])

    cells = []
    for ix in range(len(xs) - 1):
        Vx, Fx = clip_slab(V, F, 0, xs[ix], xs[ix + 1])
        if len(Fx) == 0: continue
        for iy in range(len(ys) - 1):
            Vy, Fy = clip_slab(Vx, Fx, 1, ys[iy], ys[iy + 1])
            if len(Fy) == 0: continue
            for iz in range(len(zs) - 1):
                Vc, Fc = clip_slab(Vy, Fy, 2, zs[iz], zs[iz + 1])
                if len(Fc) >= 4:
                    if do_finalize: Vc, Fc = finalize(Vc, Fc)
                    if len(Fc) >= 4:
                        cells.append((ix, iy, iz, Vc, Fc))
    return cells

def preview_pieces(cells, gap=0.5, title="크기(mm) 지정 절단 분해도"):
    if not cells:
        print("표시할 조각이 없습니다."); return
    mn = np.array([1e9]*3); mx = np.array([-1e9]*3)
    for *_, V, _ in cells:
        mn = np.minimum(mn, V.min(0)); mx = np.maximum(mx, V.max(0))
    ctr = (mn + mx) / 2; span = (mx - mn)
    fig = plt.figure(figsize=(7, 7))
    ax = fig.add_subplot(111, projection="3d")
    cmap = plt.get_cmap("tab20")
    for k, (ix, iy, iz, V, F) in enumerate(cells):
        off = (V.mean(0) - ctr) / (span/2 + 1e-9) * span * gap
        ax.add_collection3d(Poly3DCollection(
            V[F] + off, facecolor=cmap((k % 20)/20),
            edgecolor=(0, 0, 0, 0.25), linewidths=0.15))
    R = span.max() * (1 + gap) * 0.75
    ax.set_xlim(ctr[0]-R, ctr[0]+R); ax.set_ylim(ctr[1]-R, ctr[1]+R); ax.set_zlim(ctr[2]-R, ctr[2]+R)
    ax.set_box_aspect((1, 1, 1)); ax.set_axis_off()
    allwt = all(watertight_report(c[3], c[4])["watertight"] for c in cells)
    ax.set_title(f"{title}\n{len(cells)} 조각 · 전부 watertight: {allwt}")
    ax.view_init(elev=22, azim=35)
    plt.tight_layout(); plt.show()

# ----------------------------------------------------------------------------
# 3. 모델 처리 실행
# ----------------------------------------------------------------------------
if not os.path.exists(input_file_path):
    print(f"오류: 입력 파일 '{input_file_path}'을(를) 찾을 수 없습니다. 경로를 다시 확인해주세요.")
else:
    print(f"'{input_file_path}' 파일을 로드하는 중...")
    V_main, F_main = load_mesh(input_file_path)
    fmt = os.path.splitext(input_file_path)[1].lower().replace('.', '')
    
    # 모델의 실제 물리적 크기 계산
    model_size = V_main.max(0) - V_main.min(0)
    print(f"입력 모델 상태: {len(V_main)} 정점 / {len(F_main)} 삼각형")
    print(f"모델 물리적 크기(Bounding Box): X={model_size[0]:.2f}mm, Y={model_size[1]:.2f}mm, Z={model_size[2]:.2f}mm")
    
    print(f"설정된 절단 간격: X={cut_size_x_mm}mm, Y={cut_size_y_mm}mm, Z={cut_size_z_mm}mm")
    cells = grid_cut_by_size(V_main, F_main, cut_size_x_mm, cut_size_y_mm, cut_size_z_mm)
    allwt = all(watertight_report(c[3], c[4])["watertight"] for c in cells)
    
    print(f"→ 총 {len(cells)} 조각 생성 완료 · 전부 watertight: {allwt}")

    os.makedirs(output_directory_path, exist_ok=True)
    for ix, iy, iz, Vc, Fc in cells:
        save_mesh(os.path.join(output_directory_path, f"cell_{ix}_{iy}_{iz}.{fmt}"), Vc, Fc)
    print(f"저장 완료 → {output_directory_path}/")

    if show_preview:
        preview_pieces(cells)

    if download_result_zip:
        import shutil
        zip_path = shutil.make_archive(output_directory_path, "zip", output_directory_path)
        try:
            from google.colab import files
            files.download(zip_path)
            print(f"ZIP 파일 다운로드 요청됨: {zip_path}")
        except Exception:
            print(f"로컬 다운로드를 지원하지 않는 환경입니다. 생성된 ZIP 파일: {zip_path}")
