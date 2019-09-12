#include <QtWidgets>
#include "mainwindow.h"
#include "ChildWindow.h"
#include "ImageWidget.h"
#include <iostream>
#include <qinputdialog.h>

using namespace std;

MainWindow::MainWindow(QWidget *parent)
	: QMainWindow(parent)
{
	//ui.setupUi(this);

	mdi_area_ = new QMdiArea;
	mdi_area_->setHorizontalScrollBarPolicy(Qt::ScrollBarAsNeeded);
	mdi_area_->setVerticalScrollBarPolicy(Qt::ScrollBarAsNeeded);
	setCentralWidget(mdi_area_);

	window_mapper_ = new QSignalMapper(this);
	connect(window_mapper_, SIGNAL(mapped(QWidget*)), this, SLOT(SetActiveSubWindow(QWidget*)));

	CreateActions();
	CreateMenus();
	CreateToolBars();
	CreateStatusBar();

	I_ = new QImage();
	J_ = new QImage();
	K_ = new QImage();

	image_solver_ = NULL;
}

MainWindow::~MainWindow()
{
	image_solver_ = NULL;
}

void MainWindow::TestChoose()
{
	//this part maybe wrong
	ChildWindow *window = qobject_cast<ChildWindow *>(mdi_area_->activeSubWindow()->widget());
	I_ = window->imagewidget_->image();
	//image_solver_->SetAimImage(I_);

	QString filename = QFileDialog::getOpenFileName(this);

	if (!filename.isEmpty())
	{
		K_->load(filename);
	}

	if (I_->height() != K_->height() || I_->width() != K_->width()) {
		cout << "This is not the image corresponding!" << endl;
		abort();
	}

	window->imagewidget_->SetChooseforTest();
	choose_status_ = kTest;
}

void MainWindow::CreateActions()
{
	// 	action_new_ = new QAction(QIcon(":/MainWindow/Resources/images/new.png"), tr("&New"), this);
	// 	action_new_->setShortcut(QKeySequence::New);
	// 	action_new_->setStatusTip(tr("Create a new file"));

		// File IO
	action_open_ = new QAction(QIcon(":/MainWindow/Resources/images/open.png"), tr("&Open..."), this);
	action_open_->setShortcuts(QKeySequence::Open);
	action_open_->setStatusTip(tr("Open an existing file"));
	connect(action_open_, SIGNAL(triggered()), this, SLOT(Open()));

	action_save_ = new QAction(QIcon(":/MainWindow/Resources/images/save.png"), tr("&Save"), this);
	action_save_->setShortcuts(QKeySequence::Save);
	action_save_->setStatusTip(tr("Save the document to disk"));
	connect(action_save_, SIGNAL(triggered()), this, SLOT(Save()));

	action_saveas_ = new QAction(tr("Save &As..."), this);
	action_saveas_->setShortcuts(QKeySequence::SaveAs);
	action_saveas_->setStatusTip(tr("Save the document under a new name"));
	connect(action_saveas_, SIGNAL(triggered()), this, SLOT(SaveAs()));

	// Image processing
	action_invert_ = new QAction(tr("Inverse"), this);
	action_invert_->setStatusTip(tr("Invert all pixel value in the image"));
	connect(action_invert_, SIGNAL(triggered()), this, SLOT(Invert()));

	action_mirror_ = new QAction(tr("Mirror"), this);
	action_mirror_->setStatusTip(tr("Mirror image vertically or horizontally"));
	connect(action_mirror_, SIGNAL(triggered()), this, SLOT(Mirror()));

	action_gray_ = new QAction(tr("Grayscale"), this);
	action_gray_->setStatusTip(tr("Gray-scale map"));
	connect(action_gray_, SIGNAL(triggered()), this, SLOT(GrayScale()));

	action_restore_ = new QAction(tr("Restore"), this);
	action_restore_->setStatusTip(tr("Show origin image"));
	connect(action_restore_, SIGNAL(triggered()), this, SLOT(Restore()));

	// Poisson image editting
	action_choose_polygon_ = new QAction(QIcon(":/MainWindow/Resources/images/poly.png"), tr("PolyChoose"), this);
	action_choose_polygon_->setStatusTip(tr("Choose Polygon Region"));
	connect(action_choose_polygon_, SIGNAL(triggered()), this, SLOT(ChoosePoly()));

	action_choose_rectangle_ = new QAction(QIcon(":/MainWindow/Resources/images/rect.png"), tr("RectChoose"), this);
	action_choose_rectangle_->setStatusTip(tr("Choose Rectangle Region"));
	connect(action_choose_rectangle_, SIGNAL(triggered()), this, SLOT(ChooseRect()));

	action_choose_freely_ = new QAction(QIcon(":/MainWindow/Resources/images/free.png"), tr("FreeChoose"), this);
	action_choose_freely_->setStatusTip(tr("Choose Region Freely"));
	connect(action_choose_freely_, SIGNAL(triggered()), this, SLOT(ChooseFree()));

	action_choose_for_test_ = new QAction(QIcon(":/MainWindow/Resources/images/debug.png"), tr("TestChoose"), this);
	action_choose_for_test_->setStatusTip(tr("Choose for Test"));
	connect(action_choose_for_test_, SIGNAL(triggered()), this, SLOT(TestChoose()));

	action_paste_poisson_ = new QAction(QIcon(":/MainWindow/Resources/images/poisson.png"), tr("Traditional Solve"), this);
	action_paste_poisson_->setStatusTip(tr("Traditionnal Way to Solve Poisson Image Editting"));
	connect(action_paste_poisson_, SIGNAL(triggered()), this, SLOT(PastePoisson()));

	action_paste_squarepieces_ = new QAction(QIcon(":/MainWindow/Resources/images/square.png"), tr("SquarePieces Solve"), this);
	action_paste_squarepieces_->setStatusTip(tr("Solve Poisson Image Editting in Our Way"));
	connect(action_paste_squarepieces_, SIGNAL(triggered()), this, SLOT(PasteSquarePiece()));

	action_start_solving_ = new QAction(QIcon(":/MainWindow/Resources/images/start.png"), tr("Start Solving"), this);
	action_start_solving_->setStatusTip(tr("Start Solving"));
	connect(action_start_solving_, SIGNAL(triggered()), this, SLOT(StartSolve()));
}

void MainWindow::CreateMenus()
{
	menu_file_ = menuBar()->addMenu(tr("&File"));
	menu_file_->setStatusTip(tr("File menu"));
	//	menu_file_->addAction(action_new_);
	menu_file_->addAction(action_open_);
	menu_file_->addAction(action_save_);
	menu_file_->addAction(action_saveas_);

	//	menu_file_->addAction(action_choose_polygon_);

	menu_edit_ = menuBar()->addMenu(tr("&Edit"));
	menu_edit_->setStatusTip(tr("Edit menu"));
	menu_edit_->addAction(action_invert_);
	menu_edit_->addAction(action_mirror_);
	menu_edit_->addAction(action_gray_);
	menu_edit_->addAction(action_restore_);
}

void MainWindow::CreateToolBars()
{
	toolbar_file_ = addToolBar(tr("File"));
	//	toolbar_file_->addAction(action_new_);
	toolbar_file_->addAction(action_open_);
	toolbar_file_->addAction(action_save_);

	// Add separator in toolbar 
	toolbar_file_->addSeparator();
	toolbar_file_->addAction(action_invert_);
	toolbar_file_->addAction(action_mirror_);
	toolbar_file_->addAction(action_gray_);
	toolbar_file_->addAction(action_restore_);

	// Poisson Image Editing
	toolbar_file_->addSeparator();
	toolbar_file_->addAction(action_choose_rectangle_);
	toolbar_file_->addAction(action_choose_polygon_);
	toolbar_file_->addAction(action_choose_freely_);
	toolbar_file_->addAction(action_choose_for_test_);
	toolbar_file_->addAction(action_paste_poisson_);
	toolbar_file_->addAction(action_paste_squarepieces_);
	toolbar_file_->addAction(action_start_solving_);
}

void MainWindow::CreateStatusBar()
{
	statusBar()->showMessage(tr("Ready"));
}

void MainWindow::Open()
{
	QString filename = QFileDialog::getOpenFileName(this);
	if(!filename.isEmpty())
	{
		QMdiSubWindow *existing = FindChild(filename);

		if (existing)
		{
			mdi_area_->setActiveSubWindow(existing);
			return;
		}

		ChildWindow *child = CreateChildWindow();
		if (child->LoadFile(filename))
		{
			statusBar()->showMessage(tr("File loaded"), 2000);
			child->show();

			// Change size of the child window so it can fit image size
			int x = child->geometry().x();
			int y = child->geometry().y();
			int width = child->imagewidget_->ImageWidth();
			int height = child->imagewidget_->ImageHeight();
			mdi_area_->activeSubWindow()->setFixedSize(width+2*x, height+x+y);
		}
		else
		{
			child->close();
		}
	}
}

void MainWindow::Save()
{
	SaveAs();
}

void MainWindow::SaveAs()
{
	ChildWindow *window = qobject_cast<ChildWindow *>(mdi_area_->activeSubWindow()->widget());
	window->imagewidget_->SaveAs();
}

ChildWindow *MainWindow::CreateChildWindow()
{
	ChildWindow *child = new ChildWindow;
	mdi_area_->addSubWindow(child);

	return child;
}

void MainWindow::SetActiveSubWindow(QWidget* window)
{
	if (!window)
	{
		return;
	}

	mdi_area_->setActiveSubWindow(qobject_cast<QMdiSubWindow *>(window));
}

void MainWindow::Invert()
{
	ChildWindow *window = qobject_cast<ChildWindow *>(mdi_area_->activeSubWindow()->widget());
	window->imagewidget_->Invert();
}

void MainWindow::Mirror()
{
	ChildWindow *window = qobject_cast<ChildWindow *>(mdi_area_->activeSubWindow()->widget());
	window->imagewidget_->Mirror();
}

void MainWindow::GrayScale()
{
	ChildWindow *window = qobject_cast<ChildWindow *>(mdi_area_->activeSubWindow()->widget());
	window->imagewidget_->TurnGray();
}

void MainWindow::Restore()
{
	ChildWindow *window = qobject_cast<ChildWindow *>(mdi_area_->activeSubWindow()->widget());
	window->imagewidget_->Restore();
}

void MainWindow::ChooseRect()
{
	// Set source child window
	ChildWindow *window = qobject_cast<ChildWindow *>(mdi_area_->activeSubWindow()->widget());
	window->imagewidget_->set_draw_status_to_choose();
	window->imagewidget_->SetChooseRectangle();

	rect_ = window->imagewidget_->rectangle_;

	I_ = window->imagewidget_->image();

	choose_status_ = kRect;
}

void MainWindow::ChoosePoly()
{
	// Set source child window
	ChildWindow *window = qobject_cast<ChildWindow *>(mdi_area_->activeSubWindow()->widget());
	window->imagewidget_->set_draw_status_to_choose();
	window->imagewidget_->SetChoosePolygon();

	int w = window->imagewidget_->image()->width();
	int h = window->imagewidget_->image()->height();

	poly_ = window->imagewidget_->polygon_;

	I_ = window->imagewidget_->image();

	choose_status_ = kPoly;
}

void MainWindow::ChooseFree()
{
	// Set source child window
	ChildWindow *window = qobject_cast<ChildWindow *>(mdi_area_->activeSubWindow()->widget());
	window->imagewidget_->set_draw_status_to_choose();
	window->imagewidget_->SetChooseFree();

	poly_ = window->imagewidget_->polygon_;

	I_ = window->imagewidget_->image();

	choose_status_ = kFree;
}

void MainWindow::PastePoisson()
{
	method_status_ = mPS;
	// Paste image rect region to object image
	ChildWindow *child = qobject_cast<ChildWindow *>(mdi_area_->activeSubWindow()->widget());
	//child->imagewidget_->set_draw_status_to_paste();
	J_ = child->imagewidget_->image();

	image_solver_ = new PoissonSolver;

	child->imagewidget_->SetPaste(choose_status_);
	if (choose_status_ == kRect) {
		child->imagewidget_->rectangle_->paste(rect_);
		image_solver_->ImageReigstInnerVetrtice(I_, *rect_);
	}

	if (choose_status_ == kPoly || choose_status_ == kFree) {
		child->imagewidget_->polygon_->paste(poly_);
		image_solver_->ImageReigstInnerVetrtice(I_, *poly_);
	}

	if (choose_status_ = kTest) {
		image_solver_->ImageRegistforTest(I_, K_);
	}

	child->imagewidget_->set_source_window(child);
}

void MainWindow::PasteSquarePiece()
{
	method_status_ = mSP;
	// Paste image rect region to object image
	ChildWindow *child = qobject_cast<ChildWindow *>(mdi_area_->activeSubWindow()->widget());
	J_ = child->imagewidget_->image();
	child->imagewidget_->SetPaste(choose_status_);
	image_solver_ = new SquarePieces;
	if (choose_status_ == kRect) {
		child->imagewidget_->rectangle_->paste(rect_);
		image_solver_->ImageReigstInnerVetrtice(I_, *rect_);
		statusBar()->showMessage("Regist of interior is all right.");
	}
	if (choose_status_ == kPoly || choose_status_ == kFree) {
		child->imagewidget_->polygon_->paste(poly_);
		image_solver_->ImageReigstInnerVetrtice(I_, *poly_);
		statusBar()->showMessage("Regist of interior is all right.");
	}
	if (choose_status_ == kTest) {
		image_solver_->ImageRegistforTest(I_, K_);
		statusBar()->showMessage("Regist of interior is all right.");
	}

	child->imagewidget_->set_source_window(child);
}

void MainWindow::StartSolve()
{
	ChildWindow *child = qobject_cast<ChildWindow *>(mdi_area_->activeSubWindow()->widget());

	child->imagewidget_->SetDrawDefault();


	if (method_status_ == mPS && choose_status_ != kTest) {
		image_solver_->ImageRegistBoundary(J_, child->imagewidget_->point_transform_.x(), child->imagewidget_->point_transform_.y());
		statusBar()->showMessage("Regist of boundary is all right.");
		image_solver_->ParameterSolve();
		statusBar()->showMessage("The solve of Parameter is all right.");
		image_solver_->Reconstruction();
		statusBar()->showMessage("Reconstruction is all right.");

	}
	if (method_status_ == mSP && choose_status_ != kTest) {
		image_solver_->ImageRegistBoundary(J_, child->imagewidget_->point_transform_.x(), child->imagewidget_->point_transform_.y());
		statusBar()->showMessage("Regist of boundary is all right.");
		image_solver_->ParameterSolve();
		statusBar()->showMessage("The solve of Parameter is all right.");
		image_solver_->Reconstruction();
		statusBar()->showMessage("Reconstruction is all right.");
	}
	if (choose_status_ == kTest) {
		int jx, jy;
		jx = QInputDialog::getInt(this, tr("J's x coordinate"), tr("Input J's x coordinate"), 63, 0, J_->width());
		jy = QInputDialog::getInt(this, tr("J's y coordinate"), tr("Input J's y coordinate"), 195, 0, J_->height());
		image_solver_->ImageRegistBoundary(J_, jx, jy);
		statusBar()->showMessage("Regist of boundary is all right.");
		image_solver_->ParameterSolve();
		statusBar()->showMessage("The solve of Parameter is all right.");
		image_solver_->Reconstruction();
		statusBar()->showMessage("Reconstruction is all right.");
		image_solver_->SaveData();
	}

	child->imagewidget_->SetChooseDefault();
	child->imagewidget_->update();
}

void MainWindow::Paste()
{
	// Paste image rect region to object image
	ChildWindow *child = qobject_cast<ChildWindow *>(mdi_area_->activeSubWindow()->widget());
	child->imagewidget_->set_draw_status_to_paste();
	child->imagewidget_->set_source_window(child_source_);
}

QMdiSubWindow *MainWindow::FindChild(const QString &filename)
{
	QString canonical_filepath = QFileInfo(filename).canonicalFilePath();

	foreach (QMdiSubWindow *window, mdi_area_->subWindowList())
	{
		ChildWindow *child = qobject_cast<ChildWindow *>(window->widget());
		if (child->current_file() == canonical_filepath)
		{
			return window;
		}
	}

	return 0;
}
