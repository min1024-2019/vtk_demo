// vtkPickPoint.cpp : Defines the entry point for the console application.
//
#include "stdafx.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkProperty.h"
#include "vtkCubeSource.h"
#include "vtkPolyDataMapper.h"
#include "vtkActor.h"
#include "vtkCamera.h"
#include "vtkLight.h"
#include "vtkSmartPointer.h"
#include "vtkPLYReader.h"
#include "vtkPolyData.h"
#include "vtkLight.h"
#include "vtkKdTreePointLocator.h"
#include "vtkPointdata.h"
#include "vtkPolyLine.h"
#include "vtkCellArray.h"
#include "vtkCoordinate.h"
#include "vtkPolyDataMapper2D.h"
#include "vtkActor2D.h"
#include "vtkProperty2D.h" 
#include "vtkParametricSpline.h"
#include "vtkParametricFunctionSource.h"
#include "vtkTexture.h"
#include "vtkJPEGReader.h"

// for pick
#include "vtkCellPicker.h"
#include "vtkInteractorStyleTrackballCamera.h"
#include "vtkCallbackCommand.h"
#include "vtkTextMapper.h"
#include "vtkActor2D.h"
#include "vtkTextProperty.h"

#include "vtkSphereSource.h"
#include "vtkLODActor.h"
#include "vtkConeSource.h"
#include "vtkGlyph3D.h"

// for draw point
#include "vtkPoints.h"
#include "vtkVertexGlyphFilter.h"

#include "vtkTubeFilter.h"
#include "vtkDoubleArray.h"

#include <math.h>
#include <ctime>

#include <string>
#include <vector>

//const int POINTSCOUNT = 15;
const float CLOSED_DISTANCE = 1.6;
const int TIME_GAP_ON_DOUBLECLK = 250;
const int MOSUE_MOTION_LIMIT = 2;
using std::string;	

typedef struct Point
{
	double x;
	double y;
	double z;
} MYPOINT;

MYPOINT pickedPt ;
std::vector<MYPOINT> slectPts;

void drawPoints(double x, double y, double z);
void DrawPolygon();

// pick cell
int MouseMotion;
int BePickEnd;
bool BeStarted;
bool BePreview;
bool BeRotate;
bool BeEdit;

clock_t T1;
vtkRenderer *render;
vtkRenderWindow *renWin;
//vtkWin32OpenGLRenderWindow *renWinOpengl;
vtkRenderWindowInteractor *iren;
vtkCellPicker *picker;
vtkActor2D *textActor;
vtkTextMapper *textMapper; 
// draw pt
vtkActor *ptActor;
vtkPolyDataMapper* mapperPt;
vtkPoints *points;
vtkPolyData *pointsPolydata;
vtkVertexGlyphFilter *vertexFilter;
vtkPolyData *polydataPts;

// draw polyline
vtkPolyDataMapper *splinemapper;
vtkActor *splineActor;
vtkPoints *splinePts;
vtkParametricSpline *spline; 
vtkParametricFunctionSource *functionSource;
vtkKdTreePointLocator *kdTreePointLoc;

void DrawCube();
void ReadPLY(string& strFileName);
void PickPointer(string& filename);

class PickCommand : public vtkCommand
{
public:

    static PickCommand *New() { return new PickCommand; }
    void Delete() { delete this; }
 
    virtual void Execute(vtkObject *caller, unsigned long l, void *callData)
    {
		 if (!BePickEnd && !BeEdit) // pick point, left click
    {
        if (picker->GetCellId() < 0 )
        {
            //textActor->VisibilityOff();
        }
        else
        {
            double pickPos[3];
            picker->GetPickPosition( pickPos );

            pickedPt.x = pickPos[0];
            pickedPt.y = pickPos[1];
            pickedPt.z = pickPos[2];
            
            if(!BePreview)
            {
                slectPts.push_back(pickedPt);
            }

            // for real time draw pt update
            mapperPt->RemoveInputConnection(0,polydataPts->GetProducerPort());
            polydataPts->ReleaseData();
            vertexFilter->RemoveInputConnection(0,pointsPolydata->GetProducerPort());
            pointsPolydata->ReleaseData();
            points->Reset();

            // for real time draw spline update
            splinemapper->RemoveInputConnection(0, functionSource->GetOutputPort());
            if(functionSource != NULL)
            {
                functionSource->Delete();
                functionSource = vtkParametricFunctionSource::New();
            }

            if(spline != NULL)
            {
                spline->Delete();
                spline = vtkParametricSpline::New();
            }

            splinePts->Reset();

            for(std::vector<MYPOINT>::iterator it = slectPts.begin() ; it != slectPts.end(); it++)
            {
                MYPOINT pt = (*it);
                points->InsertNextPoint(pt.x, pt.y, pt.z); // pt
                splinePts->InsertNextPoint(pt.x, pt.y, pt.z); // spline
            }

            if (BePreview)
            {
                splinePts->InsertNextPoint(pickedPt.x,pickedPt.y,pickedPt.z);
            }
        }
    }
    else if(BePickEnd && !BeEdit) // last point, double click
    {
        MYPOINT pt = slectPts.back();
        double testPoint[3] = {pt.x,pt.y,pt.z};
        MYPOINT pt0 = slectPts.front();

        double distantSqure = sqrt((pt.x - pt0.x)*(pt.x - pt0.x) + (pt.y - pt0.y)*(pt.y - pt0.y) + (pt.z - pt0.z)*(pt.z - pt0.z));
        if (distantSqure < CLOSED_DISTANCE)
        {
            slectPts.pop_back();
            slectPts.push_back(pt0);
        }

        // reset value
        mapperPt->RemoveInputConnection(0,polydataPts->GetProducerPort());
        polydataPts->ReleaseData();
        vertexFilter->RemoveInputConnection(0,pointsPolydata->GetProducerPort());
        pointsPolydata->ReleaseData();
        points->Reset();

        // for real time draw spline update
        splinemapper->RemoveInputConnection(0,functionSource->GetOutputPort());
        if(functionSource != NULL)
        {
            functionSource->Delete();
            functionSource = vtkParametricFunctionSource::New();
        }

        if(spline != NULL)
        {
            spline->Delete();
            spline = vtkParametricSpline::New();
        }

        splinePts->Reset();
        
        for(std::vector<MYPOINT>::iterator it = slectPts.begin() ; it != slectPts.end(); it++)
        {
            MYPOINT pt = (*it);

            points->InsertNextPoint(pt.x, pt.y, pt.z);	// pt

            splinePts->InsertNextPoint(pt.x, pt.y, pt.z); // spline
        }
    }
    else if (BeEdit && BePickEnd) // edit spline, after right click, then left click
    {
        // edit spline, replace the nearest point with the new one

        // find the nearest point with this pick
        double pickPos[3];
        picker->GetPickPosition( pickPos );

        pickedPt.x = pickPos[0];
        pickedPt.y = pickPos[1];
        pickedPt.z = pickPos[2];

        kdTreePointLoc->SetDataSet(pointsPolydata);
        kdTreePointLoc->BuildLocator();

        double selectPoint[3] = {pickedPt.x, pickedPt.y, pickedPt.z};

        vtkIdType iD = kdTreePointLoc->FindClosestPoint(selectPoint);
        int posit = (int) iD;
        MYPOINT ptclosed = slectPts[iD];
        //for()
        slectPts.erase(slectPts.begin()+posit);
        slectPts.insert(slectPts.begin()+posit, pickedPt);

        // reset value
        mapperPt->RemoveInputConnection(0,polydataPts->GetProducerPort());
        polydataPts->ReleaseData();
        vertexFilter->RemoveInputConnection(0,pointsPolydata->GetProducerPort());
        pointsPolydata->ReleaseData();
        points->Reset();

        // for real time draw spline update
        splinemapper->RemoveInputConnection(0,functionSource->GetOutputPort());
        if(functionSource != NULL)
        {
            functionSource->Delete();
            functionSource = vtkParametricFunctionSource::New();
        }

        if(spline != NULL)
        {
            spline->Delete();
            spline = vtkParametricSpline::New();
        }

        splinePts->Reset();

        for(std::vector<MYPOINT>::iterator it = slectPts.begin() ; it != slectPts.end(); it++)
        {
            MYPOINT pt = (*it);

            points->InsertNextPoint(pt.x, pt.y, pt.z); // pt

            splinePts->InsertNextPoint(pt.x, pt.y, pt.z); // spline
        }

    }

    pointsPolydata->SetPoints(points);
    vertexFilter->SetInputConnection(pointsPolydata->GetProducerPort());
    vertexFilter->Update();
    polydataPts->ShallowCopy(vertexFilter->GetOutput());
    mapperPt->SetInputConnection(polydataPts->GetProducerPort());

    // Visualization
	ptActor->GetProperty()->SetPointSize(8);
    ptActor->GetProperty()->LightingOn();
    ptActor->GetProperty()->SetColor(1.0,0.0,0.0);
    ptActor->VisibilityOn();

    if(splinePts->GetNumberOfPoints() >=2)
    {
        // draw spline
        spline->SetPoints(splinePts);
        //if (BePickEnd)
        //{
        //    m_spline->SetClosed(1);
        //}
        functionSource->SetParametricFunction(spline);
		//functionSource->SetUResolution(100);
        functionSource->Update();

        //splinemapper->SetInputConnection(functionSource->GetOutputPort());

		// modify for tube
		vtkDoubleArray* tubeRadius = vtkDoubleArray::New();
		unsigned int n = functionSource->GetOutput()->GetNumberOfPoints();
		 tubeRadius->SetNumberOfTuples(n);
		 tubeRadius->SetName("TubeRadius");

		 vtkPolyData* tubePolyData = vtkPolyData::New();
		 tubePolyData = functionSource->GetOutput();
		 tubePolyData->GetPointData()->AddArray(tubeRadius);
		 tubePolyData->GetPointData()->SetActiveScalars("TubeRadius");

		 // set color
		 vtkUnsignedCharArray* colors = vtkUnsignedCharArray::New();
		 colors->SetName("Colors");
		 colors->SetNumberOfComponents(3);
		 colors->SetNumberOfTuples(n);
		 for (int i = 0; i < n ;i++)
		 {
			 colors->InsertTuple3(i, 0,254,0);
		 }
		 tubePolyData->GetPointData()->AddArray(colors);

		// Create the tubes
		vtkTubeFilter* tuber = vtkTubeFilter::New();
		tuber->SetInputConnection(tubePolyData->GetProducerPort());
		tuber->SetRadius(.15);
		//tuber->SetNumberOfSides(200);

		//vtkSmartPointer<vtkPolyDataMapper> tubeMapper =
		// vtkSmartPointer<vtkPolyDataMapper>::New();
		//tubeMapper->SetInputConnection(tuber->GetOutputPort());
		//tubeMapper->SetScalarRange(tubePolyData->GetScalarRange());

		splinemapper->SetInputConnection(tuber->GetOutputPort());
		splinemapper->SetScalarRange(tubePolyData->GetScalarRange());
		splinemapper->ScalarVisibilityOn();
		splinemapper->SetScalarModeToUsePointFieldData();
		splinemapper->SelectColorArray("Colors");

       /* splineActor->GetProperty()->SetColor(0.0,0.0,1.0);
        splineActor->GetProperty()->SetLineWidth(1.0);
        splineActor->GetProperty()->LightingOn();*/
        splineActor->GetProperty()->SetOpacity(1.0);
        splineActor->VisibilityOff();
    }

	renWin->Render();
    }
};

void PickerInteractionCallback( vtkObject* vtkNotUsed(object),
                                       unsigned long event,
                                       void* clientdata,
                                       void* vtkNotUsed(calldata) )
{
    vtkInteractorStyleTrackballCamera * style = 
(vtkInteractorStyleTrackballCamera*)clientdata;
	 clock_t t2;
    switch( event )
    {
	case vtkCommand::LeftButtonPressEvent:
		{
			MouseMotion = 0;
			BeRotate = true;

			int *pick = iren->GetEventPosition();
			pickedPt.x = (double)pick[0];
			pickedPt.y = (double)pick[1];
			pickedPt.z = 0.0;
			slectPts.push_back(pickedPt);
			//DrawPolygon();

			style->OnLeftButtonDown();
		}
        break;
    case vtkCommand::LeftButtonReleaseEvent:
        
         t2 = clock();

         if (t2 - T1 < TIME_GAP_ON_DOUBLECLK && MouseMotion < MOSUE_MOTION_LIMIT)
         {
             // double click, end picking
            BePickEnd = true;
             int *pick = iren->GetEventPosition();
            picker->Pick((double)pick[0], (double)pick[1], 0.0, render);
         }
         else if ((MouseMotion == 0) && ((!BePickEnd)||BeEdit))
         {
            BeStarted = true;    
            BePreview = false;
            int *pick = iren->GetEventPosition();
            double x = (double)pick[0];
            double y = (double)pick[1];
            picker->Pick((double)pick[0], (double)pick[1], 0.0, render);
         }
         T1 = t2;
      
        BeRotate = false;
        style->OnLeftButtonUp();
        break;
    case vtkCommand::RightButtonPressEvent:
        if (BePickEnd)
        {
             BeEdit = true;
        }
       
		style->OnMiddleButtonDown();
        break;

	case vtkCommand::RightButtonReleaseEvent:
		{
			int j = 0;
			style->OnMiddleButtonUp();
		}
		break;
    case vtkCommand::MouseMoveEvent:
        ++MouseMotion;
        if (BeStarted && !BePickEnd && !BeRotate)
        {
            BePreview = true;
            int *pick = iren->GetEventPosition();
            picker->Pick((double)pick[0], (double)pick[1], 0.0, render);
        }
        style->OnMouseMove();
        break;
  /*  case vtkCommand::LeftButtonPressEvent:
        MouseMotion = 0;
        style->OnLeftButtonDown();
        break;
    case vtkCommand::LeftButtonReleaseEvent:
        if (MouseMotion == 0 && BePickEnd < 1)
        {
            int *pick = iren->GetEventPosition();
            picker->Pick((double)pick[0], (double)pick[1], 0.0, render);
        }
        style->OnLeftButtonUp();
        break;
	case vtkCommand::RightButtonPressEvent:
		BePickEnd = 1;
		style->OnRightButtonDown();
		break;
	case vtkCommand::RightButtonReleaseEvent:
		break;
    case vtkCommand::MouseMoveEvent:
        MouseMotion = 1;
        style->OnMouseMove();
        break;*/
    }
}



int _tmain(int argc, _TCHAR* argv[])
{ 	 
	// 1 draw basic graphics
	//DrawCube(); 
   
	string filename("D:\\VTK Introduction\\demo\\vtk_projects\\demo_file\\mesh.ply");	
	
	// 2 read tooth ply file
	//ReadPLY(filename);

	// 3 pick point on mesh
	PickPointer(filename);

	return 0;
}

void DrawPolygon()
{
	if( slectPts.size() == 2)
	{
		// draw line
		vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
		/*	for(unsigned int j = 0; j < slectPts.size(); j++)
		{
			points->InsertNextPoint(slectPts.at(j).x, slectPts.at(j).y, 0.0);
		}*/

		points->InsertNextPoint(2.0, 100.0, 2.0);
		points->InsertNextPoint(200.0, 12.0, 20.0);		

		vtkSmartPointer<vtkPolyLine> polyLine = vtkSmartPointer<vtkPolyLine>::New();
		polyLine->GetPointIds()->SetNumberOfIds(2);
		for(unsigned int i = 0; i < slectPts.size(); i++)
		{
			polyLine->GetPointIds()->SetId(i,i);
		}

		// Create a cell array to store the lines in and add the lines to it
		vtkSmartPointer<vtkCellArray> cells = vtkSmartPointer<vtkCellArray>::New();
		cells->InsertNextCell(polyLine);

		// Create a polydata to store everything in
		vtkSmartPointer<vtkPolyData> polyData = vtkSmartPointer<vtkPolyData>::New();

		// Add the points to the dataset
		polyData->SetPoints(points);

		// Add the lines to the dataset
		polyData->SetLines(cells);

		//splinemapper = vtkPolyDataMapper::New();
		//splineActor = vtkActor::New();
		//splineActor->VisibilityOff();
		//splineActor->SetMapper(splinemapper);

		// Setup actor and mapper
		vtkSmartPointer<vtkPolyDataMapper2D> mapper = vtkSmartPointer<vtkPolyDataMapper2D>::New();
		//vtkSmartPointer<vtkCoordinate> coor = vtkSmartPointer<vtkCoordinate>::New();
		vtkSmartPointer<vtkActor2D> actor = vtkSmartPointer<vtkActor2D>::New();
		//coor->SetCoordinateSystemToDisplay();
		
		mapper->SetInput(polyData);
		//mapper->SetTransformCoordinate(coor);
		actor->SetMapper(mapper);
		actor->GetProperty()->SetColor(1.0, 0.0, 0.0);
		actor->GetProperty()->SetLineWidth(6.0);
		
		actor->VisibilityOn();
		//actor->GetProperty()->SetColor(1.0, 0.0, 0.0);
		render->AddActor(actor);
		render->Render();
		//render->AddActor(splineActor);
		//render->WorldToDisplay();	

	}


}

void DrawCube()
{
	// create render
	vtkRenderWindow *renWin = vtkRenderWindow::New();
	vtkRenderer *ren = vtkRenderer::New();
	renWin->AddRenderer(ren);

	// create iren
	vtkRenderWindowInteractor *iren = vtkRenderWindowInteractor::New();
	iren->SetRenderWindow(renWin);	

	// create cube data and mapper
	vtkCubeSource *Cube = vtkCubeSource::New();
	vtkPolyDataMapper *cubeMapper = vtkPolyDataMapper::New();
	cubeMapper->SetInput(Cube->GetOutput());

	// set actor
	vtkProperty *prop = vtkProperty::New();
	prop->SetColor(0.9804, 0.5020, 0.4471);
	prop->SetAmbient(1.0);
	prop->SetDiffuse(0.0);

	// set cube
	vtkActor *cubeActor = vtkActor::New();
	cubeActor->SetMapper(cubeMapper);
	cubeActor->SetProperty(prop);
	cubeActor->SetPosition(0, 0, 0);
	cubeActor->SetScale(1, 1, 6);

	// create camera
	vtkCamera *camera = vtkCamera::New();
	camera->SetPosition(20, 20, 50);
	camera->SetFocalPoint(0, 0, 0);
	camera->ComputeViewPlaneNormal();

	// create light
	vtkLight *light = vtkLight::New();
	light->SetPosition(50, 50, 50);
	light->SetFocalPoint(0, 0, 0);

	// create actor
	ren->AddActor(cubeActor);
	ren->SetActiveCamera(camera);
	ren->AddLight(light);
	ren->SetBackground(0, 0, 0);

	// set wind size
	renWin->SetSize(800, 600);

	// rend
	renWin->Render();
	iren->Start();

	// delete
	ren->Delete();
	renWin->Delete();
	iren->Delete();
	prop->Delete();
	Cube->Delete();
	cubeMapper->Delete();
	cubeActor->Delete();
	camera->Delete();
}

void ReadPLY(string& strFileName)
{
	// Read and display for verification
	vtkSmartPointer<vtkPLYReader> reader = vtkSmartPointer<vtkPLYReader>::New();
	reader->SetFileName(strFileName.c_str());
	reader->Update();

	vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	mapper->SetInputConnection(reader->GetOutputPort());

	vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
	actor->SetMapper(mapper);

	actor->GetProperty()->SetAmbient(1.0);
	actor->GetProperty()->SetDiffuse(0.0); 

	vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
	vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
	renderWindow->AddRenderer(renderer);
	vtkSmartPointer<vtkRenderWindowInteractor> iren = vtkSmartPointer<vtkRenderWindowInteractor>::New();
	iren->SetRenderWindow(renderWindow);

	renderer->AddActor(actor);
	renderer->SetBackground(.3, .6, .3);
	renderWindow->SetSize(1280,768);
	renderWindow->Render();
	iren->Start(); 
}

void PickPointer(string& filename){	
	// pick cell
	MouseMotion = 0;
	BePickEnd = 0;

	BeStarted = false;
	BePreview = false;
	BeRotate = false;
	BeEdit = false;
	T1 = clock();

	// read model from ply
	vtkPLYReader* reader = vtkPLYReader::New();
	reader->SetFileName(filename.c_str());
	reader->Update();

	vtkPolyDataMapper* mapper = vtkPolyDataMapper::New();
	mapper->SetInputConnection(reader->GetOutputPort());

	vtkActor* actorPly = vtkActor::New();
	actorPly->SetMapper(mapper);

	// Create a cell picker.
	PickCommand* pickObserver = PickCommand::New();
	picker = vtkCellPicker::New();
	picker->AddObserver(vtkCommand::EndPickEvent, pickObserver);

	// Create a text mapper and actor to display the results of picking.
	textMapper = vtkTextMapper::New();
	vtkTextProperty *tprop = textMapper->GetTextProperty();
	tprop->SetFontFamilyToArial();
	tprop->SetFontSize(12);
	tprop->BoldOn();
	// tprop->ShadowOn();
	tprop->SetColor(1, 0, 0);
	textActor = vtkActor2D::New();
	textActor->VisibilityOff();
	textActor->SetMapper(textMapper);

	// for draw pt
	kdTreePointLoc = vtkKdTreePointLocator::New();

	points = vtkPoints::New();
	pointsPolydata = vtkPolyData::New();
	vertexFilter = vtkVertexGlyphFilter::New();
	polydataPts = vtkPolyData::New();

	mapperPt = vtkPolyDataMapper::New();
	ptActor = vtkActor::New();
	ptActor->VisibilityOff();
	ptActor->SetMapper(mapperPt);

	//draw polyline
	spline = vtkParametricSpline::New();
	functionSource = vtkParametricFunctionSource::New();

	splinePts = vtkPoints::New();
	splinemapper = vtkPolyDataMapper::New();
	splineActor = vtkActor::New();
	splineActor->VisibilityOff();
	splineActor->SetMapper(splinemapper);

	// Create the Renderer, RenderWindow, and RenderWindowInteractor
	vtkInteractorStyleTrackballCamera *style =
		vtkInteractorStyleTrackballCamera::New();
	vtkCallbackCommand * pickerCommand = vtkCallbackCommand::New();
	pickerCommand->SetClientData(style);
	pickerCommand->SetCallback(PickerInteractionCallback);
	style->AddObserver(vtkCommand::LeftButtonPressEvent, pickerCommand);
	style->AddObserver(vtkCommand::MouseMoveEvent, pickerCommand);
	style->AddObserver(vtkCommand::LeftButtonReleaseEvent, pickerCommand);
	style->AddObserver(vtkCommand::RightButtonPressEvent, pickerCommand);
	style->AddObserver(vtkCommand::RightButtonReleaseEvent, pickerCommand);

	render = vtkRenderer::New();
	renWin = vtkRenderWindow::New();

	renWin->AddRenderer(render);
	iren = vtkRenderWindowInteractor::New();
	iren->SetRenderWindow(renWin);
	iren->SetInteractorStyle(style);
	iren->SetPicker(picker);

	// Add the actors to the renderer, set the background and size
	render->AddActor(actorPly);
	render->AddActor(ptActor);
	render->AddActor(splineActor);

	render->SetBackground(0.05, 0.05, 0.05);
	renWin->SetSize(1200, 800);
	renWin->SetCurrentCursor(10);
	iren->Initialize();
	iren->Start();

	// free
	picker->RemoveObserver(pickObserver);
	actorPly->Delete();
	picker->Delete();
	textMapper->Delete();
	textActor->Delete();
	pickerCommand->Delete();
	style->Delete();
	render->Delete();
	renWin->Delete();
	pickObserver->Delete();

	// spline free
	splinemapper->Delete();
	splineActor->Delete();
	splinePts->Delete();
	spline->Delete();
	functionSource->Delete();
	slectPts.clear();
}


