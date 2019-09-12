#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QtWidgets/QMainWindow>
#include "ui_mainwindow.h"

#include "SquarePieces.h"
#include "SquarePieces2.h"
#include "FastSquarePieces2.h"
#include "PoissonSolver.h"

#include "ImageWidget.h"

class ChildWindow;
QT_BEGIN_NAMESPACE
class QAction;
class QMenu;
class QMdiArea;
class QMdiSubWindow;
class QSignalMapper;
QT_END_NAMESPACE

enum MethodStatus {
	mDefault,
	mSP,
	mPS
};

class MainWindow : public QMainWindow
{
	Q_OBJECT

public:
	MainWindow(QWidget *parent = 0);
	~MainWindow();

private slots:
	// File IO
//	void New();
	void Open();								// Open an existing file
	void Save();								// Save image to file
	void SaveAs();
	ChildWindow *CreateChildWindow();
	void SetActiveSubWindow(QWidget* window);

	// Image Processing
	void Invert();								// Invert each pixel's rgb value
	void Mirror();								// Mirror image vertically or horizontally
	void GrayScale();							// Turn image to gray-scale map
	void Restore();								// Restore image to origin

	// Poisson Image Editing
	void Paste();								// Paste rect region to object image

	// Poisson Image Editing
	void ChooseRect();							// Choose rectangle region
	void ChoosePoly();							// Choose polygon region
	void ChooseFree();							// Choose region freely
	void PastePoisson();						// Paste rect region to object image
	void PasteSquarePiece();					// Paste rect region to object image

	void StartSolve();
	void TestChoose();

private:
	void CreateActions();
	void CreateMenus();
	void CreateToolBars();
	void CreateStatusBar();

	QMdiSubWindow *FindChild(const QString &filename);
	
private:
	Ui::MainWindowClass ui;

	QMenu						*menu_file_;
	QMenu						*menu_edit_;
	QMenu						*menu_help_;
	QToolBar					*toolbar_file_;
	QToolBar					*toolbar_edit_;
//	QAction						*action_new_;
	QAction						*action_open_;
	QAction						*action_save_;
	QAction						*action_saveas_;

	QAction						*action_invert_;
	QAction						*action_mirror_;
	QAction						*action_gray_;
	QAction						*action_restore_;

	QAction						*action_choose_rectangle_;
	QAction						*action_choose_polygon_;
	QAction						*action_choose_freely_;
	QAction						*action_choose_for_test_;

	//QAction					*action_copy_;
	QAction						*action_paste_poisson_;
	QAction						*action_paste_squarepieces_;
	QAction						*action_start_solving_;

	QMdiArea					*mdi_area_;
	QSignalMapper				*window_mapper_;

	ChildWindow					*child_source_;

	ImageSolver					*image_solver_;

	ChooseStatus				choose_status_;
	MethodStatus				method_status_;

	MyPolygon					*poly_;
	MyRectangle					*rect_;

	QImage						*I_;
	QImage						*J_;
	QImage						*K_;
};

#endif // MAINWINDOW_H
