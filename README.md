# 2D_comp_viscos
steady 2D compressible RANS Navier-Stokes solver (作成中)

つくるもの(予定)：
1.境界適合格子生成プログラム(python + fortran) → 2次元3角形格子生成プログラムに変更(2018/11/01)
2.ソルバー本体(fortran)
3.後処理(python or fortran)

※1:ソルバーはシングルスレッドで実行する前提なのでOpenMPおよびMPI並列の予定は無し → シングルスレッドを複数並べる形でジョブ自体を並列化する可能性あり
※2:一様流中におかれた物体の空力解析用

ソルバーの詳細：
離散化：セル中心型FVM
非粘性流束：SLAU2 or Roe's FDS
粘性流束：2次精度中心差分
乱流モデル：Baldwin-Lomax (for unstructured grid)
時間積分：LU-SGS, Explicit Euler (for validation)
可視化：unstructured vtk
