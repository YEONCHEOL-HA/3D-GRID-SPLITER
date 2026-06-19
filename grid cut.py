import os, struct, glob, zipfile
from collections import defaultdict
import numpy as np
from scipy.spatial import cKDTree

# ------------------------------------------------------------------ #@title 설정 (Form)
input_file_path      = "/content/nanotyrannus_scale.stl"   #@param {type:"string"}
output_directory_path= "/content/pieces1"             #@param {type:"string"}
cut_size_x_mm        = 200.0   #@param {type:"number"}
cut_size_y_mm        = 200.0   #@param {type:"number"}
cut_size_z_mm        = 200.0   #@param {type:"number"}
repair_source_holes  = True    #@param {type:"boolean"}
output_format        = "stl"   #@param ["stl", "obj"]
download_result_zip  = True    #@param {type:"boolean"}
show_preview         = True    #@param {type:"boolean"}
# ------------------------------------------------------------------

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
    tree = cKDTree(V); pairs = tree.query_pairs(r=tol, output_type='ndarray')
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
def grid_cut_by_size(V, F, sx, sy, sz):
    mn = V.min(0); mx = V.max(0); eps = (mx - mn) * 1e-6 + 1e-9
    diag = np.linalg.norm(mx - mn) or 1.0; tol = 1e-7 * diag
    def planes(s, e, step, ep):
        arr = np.arange(s - ep, e + ep, step)
        if arr[-1] < e + ep: arr = np.append(arr, e + ep)
        return arr
    xs = planes(mn[0], mx[0], sx, eps[0])
    ys = planes(mn[1], mx[1], sy, eps[1])
    zs = planes(mn[2], mx[2], sz, eps[2])
    cells = []
    for ix in range(len(xs) - 1):
        Vx, Fx = slab(V, F, 0, xs[ix], xs[ix+1], tol, diag)
        if len(Fx) == 0: continue
        for iy in range(len(ys) - 1):
            Vy, Fy = slab(Vx, Fx, 1, ys[iy], ys[iy+1], tol, diag)
            if len(Fy) == 0: continue
            for iz in range(len(zs) - 1):
                Vc, Fc = slab(Vy, Fy, 2, zs[iz], zs[iz+1], tol, diag)
                if len(Fc) >= 4: cells.append((ix, iy, iz, Vc, Fc))
    return cells

# =============================== 미리보기 ===============================
def preview_pieces(cells):
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    fig = plt.figure(figsize=(9, 7)); ax = fig.add_subplot(111, projection='3d')
    cmap = plt.get_cmap('tab20')
    for k, (ix, iy, iz, V, F) in enumerate(cells):
        step = max(1, len(F) // 4000)          # 너무 무거우면 일부만 표시
        polys = V[F[::step]]
        col = Poly3DCollection(polys, alpha=0.55, facecolor=cmap(k % 20), edgecolor='none')
        ax.add_collection3d(col)
    allV = np.vstack([c[3] for c in cells])
    mn = allV.min(0); mx = allV.max(0); ctr = (mn + mx) / 2; rng = (mx - mn).max() / 2
    ax.set_xlim(ctr[0]-rng, ctr[0]+rng); ax.set_ylim(ctr[1]-rng, ctr[1]+rng); ax.set_zlim(ctr[2]-rng, ctr[2]+rng)
    ax.set_xlabel('X'); ax.set_ylabel('Y'); ax.set_zlabel('Z')
    ax.set_title(f'{len(cells)} pieces')
    plt.tight_layout(); plt.show()

# =============================== 메인 ===============================
def main():
    assert os.path.exists(input_file_path), f"입력 파일을 찾을 수 없습니다: {input_file_path}"
    os.makedirs(output_directory_path, exist_ok=True)

    print(f"불러오는 중: {input_file_path}")
    V, F = load_mesh(input_file_path)
    rep = watertight_report(V, F)
    print(f"  정점 {len(V):,} / 면 {len(F):,} | 원본 경계모서리 {rep['boundary_edges']} "
          f"(0이 아니면 원본에 구멍이 있다는 뜻)")

    if repair_source_holes and rep['boundary_edges'] > 0:
        print("원본 구멍 메우는 중 (repair_source_holes=True)...")
        V, F = fill_holes(V, F)
        r2 = watertight_report(V, F)
        print(f"  → 경계모서리 {r2['boundary_edges']}, 비매니폴드 {r2['nonmanifold']}")

    print(f"절단 중: {cut_size_x_mm} × {cut_size_y_mm} × {cut_size_z_mm} mm ...")
    cells = grid_cut_by_size(V, F, cut_size_x_mm, cut_size_y_mm, cut_size_z_mm)
    print(f"  조각 {len(cells)}개 생성")

    saved = []
    for (ix, iy, iz, Vc, Fc) in cells:
        name = f"piece_{ix}_{iy}_{iz}.{output_format}"
        path = os.path.join(output_directory_path, name)
        save_mesh(path, Vc, Fc); saved.append(path)
        r = watertight_report(Vc, Fc)
        flag = "OK" if r['watertight'] else f"경계{r['boundary_edges']}"
        print(f"    {name:24s} 면 {len(Fc):6d}  [{flag}]")

    zip_path = None
    if download_result_zip and saved:
        zip_path = output_directory_path.rstrip("/") + "_pieces.zip"
        with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as z:
            for p in saved: z.write(p, os.path.basename(p))
        print(f"ZIP 생성: {zip_path}")
        try:
            from google.colab import files
            files.download(zip_path)
        except Exception:
            pass   # Colab 환경이 아니면 건너뜀

    if show_preview and cells:
        try: preview_pieces(cells)
        except Exception as ex: print("미리보기 생략:", ex)

    return cells, saved, zip_path

if __name__ == "__main__":
    main()
