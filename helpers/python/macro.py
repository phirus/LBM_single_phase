try: paraview.simple
except: from paraview.simple import *

paraview.simple._DisableFirstRenderCameraReset()

output_ = XMLUnstructuredGridReader( FileName=['D:\\Daten\\3D\\5\\output_0.vtu', 'D:\\Daten\\3D\\5\\output_1000.vtu', 'D:\\Daten\\3D\\5\\output_2000.vtu', 'D:\\Daten\\3D\\5\\output_3000.vtu', 'D:\\Daten\\3D\\5\\output_4000.vtu', 'D:\\Daten\\3D\\5\\output_5000.vtu', 'D:\\Daten\\3D\\5\\output_6000.vtu', 'D:\\Daten\\3D\\5\\output_7000.vtu', 'D:\\Daten\\3D\\5\\output_8000.vtu', 'D:\\Daten\\3D\\5\\output_9000.vtu', 'D:\\Daten\\3D\\5\\output_10000.vtu', 'D:\\Daten\\3D\\5\\output_11000.vtu', 'D:\\Daten\\3D\\5\\output_12000.vtu', 'D:\\Daten\\3D\\5\\output_13000.vtu', 'D:\\Daten\\3D\\5\\output_14000.vtu', 'D:\\Daten\\3D\\5\\output_15000.vtu', 'D:\\Daten\\3D\\5\\output_16000.vtu', 'D:\\Daten\\3D\\5\\output_17000.vtu', 'D:\\Daten\\3D\\5\\output_18000.vtu', 'D:\\Daten\\3D\\5\\output_19000.vtu', 'D:\\Daten\\3D\\5\\output_20000.vtu', 'D:\\Daten\\3D\\5\\output_21000.vtu', 'D:\\Daten\\3D\\5\\output_22000.vtu', 'D:\\Daten\\3D\\5\\output_23000.vtu', 'D:\\Daten\\3D\\5\\output_24000.vtu', 'D:\\Daten\\3D\\5\\output_25000.vtu', 'D:\\Daten\\3D\\5\\output_26000.vtu', 'D:\\Daten\\3D\\5\\output_27000.vtu', 'D:\\Daten\\3D\\5\\output_28000.vtu', 'D:\\Daten\\3D\\5\\output_29000.vtu', 'D:\\Daten\\3D\\5\\output_30000.vtu', 'D:\\Daten\\3D\\5\\output_31000.vtu', 'D:\\Daten\\3D\\5\\output_32000.vtu', 'D:\\Daten\\3D\\5\\output_33000.vtu', 'D:\\Daten\\3D\\5\\output_34000.vtu', 'D:\\Daten\\3D\\5\\output_35000.vtu', 'D:\\Daten\\3D\\5\\output_36000.vtu', 'D:\\Daten\\3D\\5\\output_37000.vtu', 'D:\\Daten\\3D\\5\\output_38000.vtu', 'D:\\Daten\\3D\\5\\output_39000.vtu', 'D:\\Daten\\3D\\5\\output_40000.vtu', 'D:\\Daten\\3D\\5\\output_41000.vtu', 'D:\\Daten\\3D\\5\\output_42000.vtu', 'D:\\Daten\\3D\\5\\output_43000.vtu', 'D:\\Daten\\3D\\5\\output_44000.vtu', 'D:\\Daten\\3D\\5\\output_45000.vtu', 'D:\\Daten\\3D\\5\\output_46000.vtu', 'D:\\Daten\\3D\\5\\output_47000.vtu', 'D:\\Daten\\3D\\5\\output_48000.vtu', 'D:\\Daten\\3D\\5\\output_49000.vtu', 'D:\\Daten\\3D\\5\\output_50000.vtu', 'D:\\Daten\\3D\\5\\output_51000.vtu', 'D:\\Daten\\3D\\5\\output_52000.vtu', 'D:\\Daten\\3D\\5\\output_53000.vtu', 'D:\\Daten\\3D\\5\\output_54000.vtu', 'D:\\Daten\\3D\\5\\output_55000.vtu', 'D:\\Daten\\3D\\5\\output_56000.vtu', 'D:\\Daten\\3D\\5\\output_57000.vtu', 'D:\\Daten\\3D\\5\\output_58000.vtu', 'D:\\Daten\\3D\\5\\output_59000.vtu', 'D:\\Daten\\3D\\5\\output_60000.vtu', 'D:\\Daten\\3D\\5\\output_61000.vtu', 'D:\\Daten\\3D\\5\\output_62000.vtu', 'D:\\Daten\\3D\\5\\output_63000.vtu', 'D:\\Daten\\3D\\5\\output_64000.vtu', 'D:\\Daten\\3D\\5\\output_65000.vtu', 'D:\\Daten\\3D\\5\\output_66000.vtu', 'D:\\Daten\\3D\\5\\output_67000.vtu', 'D:\\Daten\\3D\\5\\output_68000.vtu', 'D:\\Daten\\3D\\5\\output_69000.vtu', 'D:\\Daten\\3D\\5\\output_70000.vtu', 'D:\\Daten\\3D\\5\\output_71000.vtu', 'D:\\Daten\\3D\\5\\output_72000.vtu', 'D:\\Daten\\3D\\5\\output_73000.vtu', 'D:\\Daten\\3D\\5\\output_74000.vtu', 'D:\\Daten\\3D\\5\\output_75000.vtu', 'D:\\Daten\\3D\\5\\output_76000.vtu', 'D:\\Daten\\3D\\5\\output_77000.vtu', 'D:\\Daten\\3D\\5\\output_78000.vtu', 'D:\\Daten\\3D\\5\\output_79000.vtu', 'D:\\Daten\\3D\\5\\output_80000.vtu', 'D:\\Daten\\3D\\5\\output_81000.vtu', 'D:\\Daten\\3D\\5\\output_82000.vtu', 'D:\\Daten\\3D\\5\\output_83000.vtu', 'D:\\Daten\\3D\\5\\output_84000.vtu', 'D:\\Daten\\3D\\5\\output_85000.vtu', 'D:\\Daten\\3D\\5\\output_86000.vtu', 'D:\\Daten\\3D\\5\\output_87000.vtu', 'D:\\Daten\\3D\\5\\output_88000.vtu', 'D:\\Daten\\3D\\5\\output_89000.vtu', 'D:\\Daten\\3D\\5\\output_90000.vtu', 'D:\\Daten\\3D\\5\\output_91000.vtu', 'D:\\Daten\\3D\\5\\output_92000.vtu', 'D:\\Daten\\3D\\5\\output_93000.vtu', 'D:\\Daten\\3D\\5\\output_94000.vtu', 'D:\\Daten\\3D\\5\\output_95000.vtu', 'D:\\Daten\\3D\\5\\output_96000.vtu', 'D:\\Daten\\3D\\5\\output_97000.vtu', 'D:\\Daten\\3D\\5\\output_98000.vtu', 'D:\\Daten\\3D\\5\\output_99000.vtu', 'D:\\Daten\\3D\\5\\output_100000.vtu', 'D:\\Daten\\3D\\5\\output_101000.vtu', 'D:\\Daten\\3D\\5\\output_102000.vtu', 'D:\\Daten\\3D\\5\\output_103000.vtu', 'D:\\Daten\\3D\\5\\output_104000.vtu', 'D:\\Daten\\3D\\5\\output_105000.vtu', 'D:\\Daten\\3D\\5\\output_106000.vtu', 'D:\\Daten\\3D\\5\\output_107000.vtu', 'D:\\Daten\\3D\\5\\output_108000.vtu', 'D:\\Daten\\3D\\5\\output_109000.vtu', 'D:\\Daten\\3D\\5\\output_110000.vtu', 'D:\\Daten\\3D\\5\\output_111000.vtu', 'D:\\Daten\\3D\\5\\output_112000.vtu', 'D:\\Daten\\3D\\5\\output_113000.vtu', 'D:\\Daten\\3D\\5\\output_114000.vtu', 'D:\\Daten\\3D\\5\\output_115000.vtu', 'D:\\Daten\\3D\\5\\output_116000.vtu', 'D:\\Daten\\3D\\5\\output_117000.vtu', 'D:\\Daten\\3D\\5\\output_118000.vtu', 'D:\\Daten\\3D\\5\\output_119000.vtu', 'D:\\Daten\\3D\\5\\output_120000.vtu', 'D:\\Daten\\3D\\5\\output_121000.vtu', 'D:\\Daten\\3D\\5\\output_122000.vtu', 'D:\\Daten\\3D\\5\\output_123000.vtu', 'D:\\Daten\\3D\\5\\output_124000.vtu', 'D:\\Daten\\3D\\5\\output_125000.vtu', 'D:\\Daten\\3D\\5\\output_126000.vtu', 'D:\\Daten\\3D\\5\\output_127000.vtu', 'D:\\Daten\\3D\\5\\output_128000.vtu', 'D:\\Daten\\3D\\5\\output_129000.vtu', 'D:\\Daten\\3D\\5\\output_130000.vtu', 'D:\\Daten\\3D\\5\\output_131000.vtu', 'D:\\Daten\\3D\\5\\output_132000.vtu'] )

AnimationScene2 = GetAnimationScene()
AnimationScene2.EndTime = 132.0
AnimationScene2.PlayMode = 'Snap To TimeSteps'

output_.CellArrayStatus = ['Psi']

RenderView2 = GetRenderView()
RenderView2.CenterOfRotation = [40.0, 40.0, 150.0]

DataRepresentation3 = Show()
DataRepresentation3.EdgeColor = [0.0, 0.0, 0.5000076295109483]
DataRepresentation3.ColorAttributeType = 'CELL_DATA'
DataRepresentation3.SelectionCellFieldDataArrayName = 'Psi'
DataRepresentation3.ColorArrayName = ('CELL_DATA', 'Psi')
DataRepresentation3.ScalarOpacityUnitDistance = 2.579662048328137
DataRepresentation3.ScaleFactor = 30.0

a1_Psi_PVLookupTable = GetLookupTableForArray( "Psi", 1, RGBPoints=[-1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 0.0], VectorMode='Magnitude', NanColor=[0.498039, 0.498039, 0.498039], ColorSpace='HSV', ScalarRangeInitialized=1.0 )

a1_Psi_PiecewiseFunction = CreatePiecewiseFunction( Points=[-1.0, 0.0, 0.5, 0.0, 0.26351361362021186, 0.304761916399002, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0] )

DataRepresentation3.ScalarOpacityFunction = a1_Psi_PiecewiseFunction
DataRepresentation3.LookupTable = a1_Psi_PVLookupTable

a1_Psi_PVLookupTable.ScalarOpacityFunction = a1_Psi_PiecewiseFunction

RenderView2.CameraPosition = [40.0, 40.0, 769.39875929653]
RenderView2.CameraFocalPoint = [40.0, 40.0, 150.0]
RenderView2.CameraClippingRange = [314.7047717035647, 1004.439740685978]
RenderView2.CameraParallelScale = 160.31219541881399

CellDatatoPointData1 = CellDatatoPointData()

DataRepresentation4 = Show()
DataRepresentation4.EdgeColor = [0.0, 0.0, 0.5000076295109483]
DataRepresentation4.ColorAttributeType = 'POINT_DATA'
DataRepresentation4.SelectionPointFieldDataArrayName = 'Psi'
DataRepresentation4.SelectionCellFieldDataArrayName = 'Psi'
DataRepresentation4.ScalarOpacityFunction = a1_Psi_PiecewiseFunction
DataRepresentation4.ColorArrayName = ('POINT_DATA', 'Psi')
DataRepresentation4.ScalarOpacityUnitDistance = 2.579662048328137
DataRepresentation4.LookupTable = a1_Psi_PVLookupTable
DataRepresentation4.ScaleFactor = 30.0

DataRepresentation3.Visibility = 0

Slice2 = Slice( SliceType="Plane" )

Slice2.SliceOffsetValues = [0.0]
Slice2.SliceType.Origin = [40.0, 40.0, 150.0]
Slice2.SliceType = "Plane"

# toggle the 3D widget visibility.
active_objects.source.SMProxy.InvokeEvent('UserEvent', 'ShowWidget')
RenderView2.CameraClippingRange = [253.93282872565942, 1080.9390917296791]

DataRepresentation5 = Show()
DataRepresentation5.EdgeColor = [0.0, 0.0, 0.5000076295109483]
DataRepresentation5.SelectionPointFieldDataArrayName = 'Psi'
DataRepresentation5.SelectionCellFieldDataArrayName = 'Psi'
DataRepresentation5.ColorArrayName = ('POINT_DATA', 'Psi')
DataRepresentation5.LookupTable = a1_Psi_PVLookupTable
DataRepresentation5.ScaleFactor = 30.0

DataRepresentation4.Visibility = 0

RenderView2.CameraViewUp = [0.0, 0.0, 1.0]
RenderView2.CameraPosition = [-559.8080508376324, 40.0, 150.0]
RenderView2.CameraClippingRange = [358.2502677993489, 905.3262960585099]
RenderView2.CameraFocalPoint = [40.0, 40.0, 150.0]
RenderView2.CameraParallelScale = 155.24174696260025
RenderView2.CenterOfRotation = [40.0, 40.0, 150.0]

RenderView2.CameraClippingRange = [593.8099703292561, 608.8051716001969]

WriteAnimation('D:/Daten/3D/5/5.avi', Magnification=1, Quality=2, FrameRate=8.000000)


Render()
