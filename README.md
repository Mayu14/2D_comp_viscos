# 2D_comp_viscos
steady 2D compressible RANS Navier-Stokes solver (作成中)

つくるもの(予定)：
1.境界適合格子生成プログラム(python)
2.ソルバー本体(fortran)
3.後処理(python or fortran)

※1:ソルバーはシングルスレッドで実行する前提なのでOpenMPおよびMPI並列の予定は無し
※2:一様流中におかれた物体の空力解析用

ソルバーの詳細：
離散化：セル中心型FVM
非粘性流束：SLAU2
粘性流束：2次精度中心差分
乱流モデル：Baldwin-Lomax
時間積分：LU-SGS, Explicit Euler (for validation)
可視化：unstructured vtk
