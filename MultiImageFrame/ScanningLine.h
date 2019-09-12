#pragma once
#ifndef _INC_STRUCTURE
#define _INC_STRUCTURE

#pragma once
#include <vector>
#include "Polygon.h"
#include "Rectangle.h"
#include <list>
#include <functional>
#include <algorithm>

typedef struct tagEDGE
{
	double xi;
	double dx;
	int ymax;
}EDGE;

typedef struct tagMARK
{
	int x;
	int y;
	bool is_boundary;

}MARK;



class ScanningLine
{
public:
	ScanningLine();
	ScanningLine(int y_min_, int y_max_) {
		y_min = y_min_;
		y_max = y_max_;
		/*start_line = new vector < PanelPoint* > ;*/
		/*for (int i = y_min_; i <= y_max_; i++){
			PanelPoint* p = new PanelPoint;
			p = NULL;
			start_line->push_back(p);
		}*/
		slNet = new std::vector< std::vector<EDGE> >(y_max - y_min + 1);
	}
	/*int abs(int x){ return x > 0 ? x : -x; }
	void AET(MyPolygon *poly);
	void simplyfy();*/
	~ScanningLine();
	std::vector< std::vector<EDGE> >* slNet;
	/*vector<PanelPoint*>* start_line;*/
	void InitScanLineNewEdgeTable(MyPolygon& pl);
	void UpdateAetEdgeInfo(EDGE& e)
	{
		e.xi += e.dx;
	}


	void UpdateAndResortAet(std::vector<EDGE>& aet);

	void RemoveNonActiveEdgeFromAet(std::vector<EDGE>& aet, int y);

	void InsertNetListToAet(std::vector<EDGE>& net, std::vector<EDGE>& aet);
	vector<int> Fill_inner_mark(const std::vector<EDGE>& aet, int y);
	vector<vector<MARK>> ProcessScanLineFill();
	vector<vector<MARK>> boundary_expansion(vector<vector<MARK>> scaninglines);
	void UpdateScanningLine(MyPolygon& pl);
	void UpdateScanningLine(MyRectangle& rt);
	void RestScannigLine();
public:
	int y_min, y_max;
	int x_min;
	vector<vector<MARK>>* inner_mark;
};
#endif