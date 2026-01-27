import pyvista as pv
import vtk
import numpy as np
import os
import random
import warnings
import gc

# ==========================================
# [Colab 전용 설정]
# ==========================================
pv.start_xvfb()
warnings.filterwarnings("ignore")

def clip_box_and_force_cap(mesh, bounds):
    """
    [FlashPrint 스타일 + 안전장치]
    1. vtkClipClosedSurface로 평평한 뚜껑 생성을 시도합니다.
    2. 만약 실패해서 구멍이 남으면, fill_holes로 강제로 막아버립니다.
    """
    xmin, xmax, ymin, ymax, zmin, zmax = bounds

    # --- 1단계: 정밀 절단 (평평한 뚜껑 생성 시도) ---
    planes = vtk.vtkPlaneCollection()

    # 평면 정의 (박스 안쪽을 향함)
    p1 = vtk.vtkPlane(); p1.SetOrigin(xmin, 0, 0); p1.SetNormal(1, 0, 0); planes.AddItem(p1)
    p2 = vtk.vtkPlane(); p2.SetOrigin(xmax, 0, 0); p2.SetNormal(-1, 0, 0); planes.AddItem(p2)
    p3 = vtk.vtkPlane(); p3.SetOrigin(0, ymin, 0); p3.SetNormal(0, 1, 0); planes.AddItem(p3)
    p4 = vtk.vtkPlane(); p4.SetOrigin(0, ymax, 0); p4.SetNormal(0, -1, 0); planes.AddItem(p4)
    p5 = vtk.vtkPlane(); p5.SetOrigin(0, 0, zmin); p5.SetNormal(0, 0, 1); planes.AddItem(p5)
    p6 = vtk.vtkPlane(); p6.SetOrigin(0, 0, zmax); p6.SetNormal(0, 0, -1); planes.AddItem(p6)

    clipper = vtk.vtkClipClosedSurface()
    clipper.SetInputData(mesh)
    clipper.SetClippingPlanes(planes)
    clipper.SetGenerateFaces(True) # 1차 시도: 뚜껑 생성
    clipper.SetTolerance(0.01)     # 오차 범위 넉넉하게
    clipper.Update()

    result = pv.wrap(clipper.GetOutput())

    # --- 2단계: 안전장치 (강제 메움) ---
    if result.n_cells > 0:
        result = result.clean(tolerance=0.01)
        if result.n_open_edges > 0:
            result = result.fill_holes(hole_size=1000000000)
        result = result.compute_normals(auto_orient_normals=True)

    return result

def split_mesh_with_red_grid(input_path, output_dir, grid_size_mm, visualize=True):

    # 1. 파일 로드
    if not os.path.exists(input_path):
        print(f"Error: 파일을 찾을 수 없습니다 -> {input_path}")
        return

    print(f"Loading mesh: {input_path}...")
    mesh = pv.read(input_path)

    # 2. 전처리
    print("메쉬 최적화 중...")
    mesh = mesh.clean().triangulate()

    if not mesh.is_manifold:
        print("원본의 구멍을 메우는 중...")
        mesh = mesh.fill_holes(hole_size=10000000)

    mesh = mesh.compute_normals(auto_orient_normals=True)

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # 3. Grid 설정
    xmin, xmax, ymin, ymax, zmin, zmax = mesh.bounds
    gx, gy, gz = grid_size_mm

    x_ranges = np.arange(xmin, xmax, gx)
    y_ranges = np.arange(ymin, ymax, gy)
    z_ranges = np.arange(zmin, zmax, gz)

    count = 0
    total_grids = len(x_ranges) * len(y_ranges) * len(z_ranges)
    print(f"Processing... (총 {total_grids}개 영역 작업 시작)")

    # ==========================================
    # [시각화 설정 변경]
    # ==========================================
    if visualize:
        p = pv.Plotter(notebook=True)
        # 배경: 원본 메쉬를 반투명(Ghost)하게 표시
        p.add_mesh(mesh, color='lightgrey', opacity=0.3, style='surface', label='Original')
        p.add_mesh(mesh, color='black', opacity=0.1, style='wireframe')
        print("[Visual Info] 원본 위에 빨간색 절단 라인을 표시합니다.")

    # 4. 분할 루프
    for i, x_start in enumerate(x_ranges):
        for j, y_start in enumerate(y_ranges):
            for k, z_start in enumerate(z_ranges):

                x_end = min(x_start + gx, xmax)
                y_end = min(y_start + gy, ymax)
                z_end = min(z_start + gz, zmax)
                bounds = [x_start, x_end, y_start, y_end, z_start, z_end]

                try:
                    # 1. 데이터 존재 여부 확인
                    if mesh.clip_box(bounds, invert=False).n_cells == 0:
                        continue

                    # 2. 자르고 저장 (로직 유지)
                    part = clip_box_and_force_cap(mesh, bounds)

                    if part.n_cells > 0:
                        # 파일 저장
                        part = part.clean().triangulate()
                        filename = f"part_x{i}_y{j}_z{k}.stl"
                        filepath = os.path.join(output_dir, filename)
                        part.save(filepath)
                        count += 1
                        print(f"  -> {filename} 저장 완료")

                        # ==========================================
                        # [시각화: 빨간색 박스만 추가]
                        # ==========================================
                        if visualize:
                            # 잘린 조각(part) 대신, '빨간색 그리드 박스'를 그립니다.
                            grid_box = pv.Box(bounds)
                            p.add_mesh(grid_box, color='red', style='wireframe', line_width=3, opacity=1.0)

                except Exception as e:
                    print(f"  Warning: {i},{j},{k} 에러: {e}")
                    continue

                gc.collect()

    print(f"\n✅ 최종 완료! 총 {count}개 파일 생성됨.")
    print(f"저장 위치: {output_dir}")

    # 5. 결과 시각화
    if visualize and count > 0:
        p.show_grid(color='black')
        p.show_axes()
        p.set_background('white') # 깔끔한 흰색 배경
        p.camera_position = 'iso'

        # 정적 이미지로 출력
        p.show(jupyter_backend='static', window_size=[1024, 768])

# ==========================================
# [설정 구역]
# ==========================================
INPUT_STL_PATH = "/content/drive/MyDrive/nanotyrannus lanceolatus_scaling.stl"
GRID_SIZE_MM = [210, 210, 210]
OUTPUT_FOLDER = "/content/drive/MyDrive/nanotyrannosaurus_split_red_grid"

if __name__ == "__main__":
    split_mesh_with_red_grid(INPUT_STL_PATH, OUTPUT_FOLDER, GRID_SIZE_MM, visualize=True)
