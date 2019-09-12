#include "ScanningLine.h"


ScanningLine::ScanningLine()
{
	y_min = 0;
	y_max = 0;
	inner_mark = new vector<vector<MARK>>;
}


ScanningLine::~ScanningLine()
{
	RestScannigLine();
}

bool MyEdgeXiComparator(EDGE& e1, EDGE& e2)
{
	return (e1.xi <= e2.xi);
}

bool MyIsEdgeOutOfActive(EDGE e, int y)
{
	return (e.ymax == y);
}

/*void ScanningLine::AET(MyPolygon *poly){
	y_min = poly->get_min_y();
	y_max = poly->get_max_y();
	start_line = new vector < PanelPoint* >;
	for (int i = y_min; i <= y_max; i++){
		PanelPoint* p = new PanelPoint;
		p = NULL;
		start_line->push_back(p);
	}
	PanelPoint *panel;
	PanelPoint *temp;
	double new_x;
	int k,dy,dir;
	int n = poly->vertexes->size();
	if (poly->vertexes->size() > 1){
		for (int i = 0; i < n; i++){
			if ((*poly->vertexes)[i%n].ry() != (*poly->vertexes)[(i + 1) % n].ry()){
				PanelPoint* p = new PanelPoint((int)(*poly->vertexes)[i%n].rx(), ((*poly->vertexes)[(i + 1) % n].rx() - (*poly->vertexes)[i%n].rx()) / ((*poly->vertexes)[(i + 1) % n].ry() - (*poly->vertexes)[i%n].ry()), (*poly->vertexes)[(i + 1) % n].ry());
				k = (*poly->vertexes)[i%n].ry() - y_min;
				panel = (*start_line)[k];
				if (panel == NULL){
					(*start_line)[k] = p;
				}
				else{
					while (panel->get_x() < p->get_x() && panel->get_next() != NULL){
						panel = panel->get_next();
					}
					if (panel->get_next() != NULL){
						temp = panel->get_next();
						panel->set_next(p);
						p->set_next(temp);
					}
					else{
						panel->set_next(p);
					}
				}
				dy = abs((*poly->vertexes)[i%n].ry() - (*poly->vertexes)[(i + 1) % n].ry());
				for (int j = 1; j <= dy ; j++){
					new_x = j*p->get_dx() + p->get_x();
					PanelPoint *new_p = new PanelPoint((int)new_x, p->get_dx(), p->get_another_y());
					dir = ((*poly->vertexes)[(i + 1) % n].ry() - (*poly->vertexes)[i%n].ry()) / abs((*poly->vertexes)[i%n].ry() - (*poly->vertexes)[(i + 1) % n].ry());
					k = (*poly->vertexes)[i%n].ry() - y_min + j*dir;
					panel = (*start_line)[k];
					if (panel == NULL){
						(*start_line)[k] = new_p;
					}
					else{
						while (panel->get_x() < new_p->get_x() && panel->get_next() != NULL){
							panel = panel->get_next();
						}
						if (panel->get_next() != NULL){
							temp = panel->get_next();
							panel->set_next(new_p);
							new_p->set_next(temp);
						}
						else{
							panel->set_next(new_p);
						}
					}
				}
			}
		}
		simplyfy();
		for (int i = 0; i < start_line->size(); i++){
			panel = (*start_line)[i];
			if (panel != NULL){
				while (panel->get_next() != NULL){
					temp = panel;
					panel = panel->get_next();
					panel->set_token(!temp->get_token());
				}
			}
		}
	}
}

void ScanningLine::simplyfy(){
	PanelPoint *panel;
	for (int i = 0; i < start_line->size(); i++){
		panel = (*start_line)[i];
		while (panel != NULL){
			if (panel->get_next() != NULL){
				if (panel->get_dx() == panel->get_next()->get_dx() && panel->get_x() == panel->get_next()->get_x()){
					panel->set_next(panel->get_next()->get_next());
				}
			}
			panel = panel->get_next();
		}
	}
}*/

void ScanningLine::InitScanLineNewEdgeTable(MyPolygon& pl) {
	y_max = pl.get_max_y();
	y_min = pl.get_min_y();
	x_min = pl.get_min_x();
	slNet = new std::vector < std::vector<EDGE> >;
	for (int i = 0; i < y_max - y_min + 1; i++) {
		vector<EDGE> p;
		slNet->push_back(p);
	}
	EDGE e;
	for (int i = 0; i < pl.vertexes->size(); i++)
	{
		QPoint ps = (*pl.vertexes)[i];
		QPoint pe = (*pl.vertexes)[(i + 1) % pl.vertexes->size()];
		QPoint pss = (*pl.vertexes)[(i - 1 + pl.vertexes->size()) % pl.vertexes->size()];
		QPoint pee = (*pl.vertexes)[(i + 2) % pl.vertexes->size()];
		if (pe.y() != ps.y()) //不处理水平线
		{
			e.dx = double(pe.x() - ps.x()) / double(pe.y() - ps.y());
			if (pe.y() > ps.y())
			{
				e.xi = ps.x();
				if (pee.y() >= pe.y())
					e.ymax = pe.y() - 1;
				else
					e.ymax = pe.y();

				(*slNet)[(ps.y() - y_min)].push_back(e);
			}
			else
			{
				e.xi = pe.x();
				if (pss.y() >= ps.y())
					e.ymax = ps.y() - 1;
				else
					e.ymax = ps.y();
				(*slNet)[(pe.y() - y_min)].push_back(e);
			}
		}
	}
	int x = 1;
}

vector<int> ScanningLine::Fill_inner_mark(const std::vector<EDGE>& aet, int y) {
	vector<int> inner_mark;
	for (int i = 0; i < aet.size(); i = i + 2) {
		if (i + 1 < aet.size()) {
			int a = (int)aet[i].xi;
			int b = (int)aet[i + 1].xi;
			for (int j = a + 1; j < b; j++) {
				int mark = j;
				inner_mark.push_back(mark);
			}
		}
	}
	return(inner_mark);
}

void ScanningLine::InsertNetListToAet(std::vector<EDGE>& net, std::vector<EDGE>& aet)
{
	aet.insert(aet.end(), net.begin(), net.end());
	sort(aet.begin(), aet.end(), MyEdgeXiComparator);
}

vector<vector<MARK>> ScanningLine::ProcessScanLineFill() {
	vector<vector<MARK>> inner_mark;
	std::vector<EDGE> aet;
	for (int y = y_min; y <= y_max; y++)
	{
		vector<MARK> inner_mark_y;
		InsertNetListToAet((*slNet)[y - y_min], aet);
		vector<int> temp_inner_mark = Fill_inner_mark(aet, y);
		sort(aet.begin(), aet.end(), MyEdgeXiComparator);
		for (int i = 0; i < temp_inner_mark.size(); i++) {
			MARK mark;
			mark.is_boundary = false;
			mark.x = temp_inner_mark[i];
			mark.y = y;
			inner_mark_y.push_back(mark);
		}
		inner_mark.push_back(inner_mark_y);
		//FillAetScanLine(aet, y, color);
		//删除非活动边
		RemoveNonActiveEdgeFromAet(aet, y);
		//更新活动边表中每项的xi值，并根据xi重新排序
		UpdateAndResortAet(aet);
	}
	return(inner_mark);
}

vector<vector<MARK>> ScanningLine::boundary_expansion(vector<vector<MARK>> scaninglines)
{
	vector<vector<MARK>> ex_scaninglines;
	vector<MARK> front_line;
	ex_scaninglines.push_back(front_line);

	for (int y = 0; y < scaninglines.size(); y++) {
		if (scaninglines[y].size() != 0) {
			vector<MARK> inner_line;
			if (scaninglines[y][0].x == 0) {
				scaninglines[y][0].is_boundary = true;
			}
			else {
				MARK first_mark;
				first_mark.y = scaninglines[y][0].y;
				first_mark.x = scaninglines[y][0].x - 1;
				first_mark.is_boundary = true;
				inner_line.push_back(first_mark);
			}
			inner_line.insert(inner_line.end(), scaninglines[y].begin(), scaninglines[y].end());
			//this part need make diffenece if we passing the maximum in scaninling line
			MARK last_mark;
			last_mark.y = scaninglines[y][scaninglines[y].size() - 1].y;
			last_mark.x = scaninglines[y][scaninglines[y].size() - 1].x + 1;
			last_mark.is_boundary = true;
			inner_line.push_back(last_mark);
			ex_scaninglines.push_back(inner_line);
		}
	}

	vector<MARK> behind_line;
	ex_scaninglines.push_back(behind_line);

	//expansion
	for (int y = 1; y < ex_scaninglines.size() - 1; y++) {
		for (int x = 1; x < ex_scaninglines[y].size() - 1; x++) {
			MARK mark = ex_scaninglines[y][x];
			//the first and last line don't need cheek
			if (!mark.is_boundary) {
				//upper expansion
				//v_0
				vector<MARK>::iterator it = ex_scaninglines[y - 1].begin();
				if (ex_scaninglines[y - 1].size() != 0) {
					while (it != ex_scaninglines[y - 1].end()) {
						if (it->x >= mark.x - 1) {
							break;
						}
						it++;
					}
					if (it == ex_scaninglines[y - 1].end()) {
						MARK mark_0;
						mark_0.y = mark.y - 1;
						mark_0.x = mark.x - 1;
						mark_0.is_boundary = true;
						it = ex_scaninglines[y - 1].insert(it, mark_0);
					}
					else {
						if (it->x > mark.x - 1) {
							MARK mark_0;
							mark_0.y = mark.y - 1;
							mark_0.x = mark.x - 1;
							mark_0.is_boundary = true;
							it = ex_scaninglines[y - 1].insert(it, mark_0);
						}
					}
				}
				else {
					MARK mark_0;
					mark_0.y = mark.y - 1;
					mark_0.x = mark.x - 1;
					mark_0.is_boundary = true;
					it = ex_scaninglines[y - 1].insert(it, mark_0);
				}
				//v_1
				it++;
				if (it == ex_scaninglines[y - 1].end()) {
					MARK mark_1;
					mark_1.y = mark.y - 1;
					mark_1.x = mark.x;
					mark_1.is_boundary = true;
					it = ex_scaninglines[y - 1].insert(it, mark_1);
				}
				else {
					if (it->x > mark.x) {
						MARK mark_1;
						mark_1.y = mark.y - 1;
						mark_1.x = mark.x;
						mark_1.is_boundary = true;
						it = ex_scaninglines[y - 1].insert(it, mark_1);
					}
				}
				//v_2
				it++;
				if (it == ex_scaninglines[y - 1].end()) {
					MARK mark_2;
					mark_2.y = mark.y - 1;
					mark_2.x = mark.x + 1;
					mark_2.is_boundary = true;
					it = ex_scaninglines[y - 1].insert(it, mark_2);
				}
				else {
					if (it->x > mark.x + 1) {
						MARK mark_2;
						mark_2.y = mark.y - 1;
						mark_2.x = mark.x + 1;
						mark_2.is_boundary = true;
						it = ex_scaninglines[y - 1].insert(it, mark_2);
					}
				}

				//lower expansion
				//v_5
				it = ex_scaninglines[y + 1].begin();
				if (ex_scaninglines[y + 1].size() != 0) {
					while (it != ex_scaninglines[y + 1].end()) {
						if (it->x >= mark.x - 1) {
							break;
						}
						it++;
					}
					if (it == ex_scaninglines[y + 1].end()) {
						MARK mark_5;
						mark_5.y = mark.y + 1;
						mark_5.x = mark.x - 1;
						mark_5.is_boundary = true;
						it = ex_scaninglines[y + 1].insert(it, mark_5);
					}
					else {
						if (it->x > mark.x - 1) {
							MARK mark_5;
							mark_5.y = mark.y + 1;
							mark_5.x = mark.x - 1;
							mark_5.is_boundary = true;
							it = ex_scaninglines[y + 1].insert(it, mark_5);
						}
					}
				}
				else {
					MARK mark_5;
					mark_5.y = mark.y + 1;
					mark_5.x = mark.x - 1;
					mark_5.is_boundary = true;
					it = ex_scaninglines[y + 1].insert(it, mark_5);
				}
				//v_6
				it++;
				if (it == ex_scaninglines[y + 1].end()) {
					MARK mark_6;
					mark_6.y = mark.y + 1;
					mark_6.x = mark.x;
					mark_6.is_boundary = true;
					it = ex_scaninglines[y + 1].insert(it, mark_6);
				}
				else {
					if (it->x > mark.x) {
						MARK mark_6;
						mark_6.y = mark.y + 1;
						mark_6.x = mark.x;
						mark_6.is_boundary = true;
						it = ex_scaninglines[y + 1].insert(it, mark_6);
					}
				}
				//v_7
				it++;
				if (it == ex_scaninglines[y + 1].end()) {
					MARK mark_7;
					mark_7.y = mark.y + 1;
					mark_7.x = mark.x + 1;
					mark_7.is_boundary = true;
					it = ex_scaninglines[y + 1].insert(it, mark_7);
				}
				else {
					if (it->x > mark.x + 1) {
						MARK mark_7;
						mark_7.y = mark.y + 1;
						mark_7.x = mark.x + 1;
						mark_7.is_boundary = true;
						it = ex_scaninglines[y + 1].insert(it, mark_7);
					}
				}
			}
		}
	}

	for (int y = 1; y < ex_scaninglines.size() - 1; y++) {
		vector<MARK>::iterator it = ex_scaninglines[y].begin();
		while (it != ex_scaninglines[y].end()) {
			if (it->is_boundary) {
				it++;
			}
			else {
				//v_3
				vector<MARK>::iterator fit = it;
				fit--;
				if (fit->x + 1 < it->x) {
					MARK mark_3;
					mark_3.x = it->x - 1;
					mark_3.y = it->y;
					mark_3.is_boundary = true;
					it = ex_scaninglines[y].insert(it, mark_3);
					it++;
				}
				//v_4
				vector<MARK>::iterator bit = it;
				bit++;
				if (it->x + 1 < bit->x) {
					MARK mark_4;
					mark_4.x = it->x - 1;
					mark_4.y = it->y;
					mark_4.is_boundary = true;
					it = ex_scaninglines[y].insert(bit, mark_4);
				}
				it++;
			}
		}
	}
	return(ex_scaninglines);
}

void ScanningLine::UpdateScanningLine(MyPolygon& pl) {
	InitScanLineNewEdgeTable(pl);
	(*inner_mark) = boundary_expansion(ProcessScanLineFill());
}

void ScanningLine::UpdateScanningLine(MyRectangle & rt)
{
	int width = rt.get_width();
	int height = rt.get_height();
	int x = rt.get_start_x();
	int y = rt.get_start_y();
	y_min = y;
	y_max = y + height;
	x_min = rt.get_start_x();
	for (int i = 0; i < height; i++) {
		vector<MARK> lines;
		inner_mark->push_back(lines);
	}
	MARK mark_left;
	mark_left.x = x;
	mark_left.y = y;
	mark_left.is_boundary = true;
	(*inner_mark)[0].push_back(mark_left);
	for (int j = 1; j < width - 1; j++) {
		MARK mark_inner;
		mark_inner.x = x + j;
		mark_inner.y = y;
		mark_inner.is_boundary = true;
		(*inner_mark)[0].push_back(mark_inner);
	}
	MARK mark_right;
	mark_right.x = x + width;
	mark_right.y = y;
	mark_right.is_boundary = true;
	(*inner_mark)[0].push_back(mark_right);
	for (int i = 1; i < height - 1; i++) {
		MARK mark_left;
		mark_left.x = x;
		mark_left.y = y + i;
		mark_left.is_boundary = true;
		(*inner_mark)[i].push_back(mark_left);
		for (int j = 1; j < width - 1; j++) {
			MARK mark_inner;
			mark_inner.x = x + j;
			mark_inner.y = y + i;
			mark_inner.is_boundary = false;
			(*inner_mark)[i].push_back(mark_inner);
		}
		MARK mark_right;
		mark_right.x = x + width;
		mark_right.y = y + i;
		mark_right.is_boundary = true;
		(*inner_mark)[i].push_back(mark_right);
	}
	//MARK mark_left;
	mark_left.x = x;
	mark_left.y = y + height;
	mark_left.is_boundary = true;
	(*inner_mark)[height - 1].push_back(mark_left);
	for (int j = 1; j < width - 1; j++) {
		MARK mark_inner;
		mark_inner.x = x + j;
		mark_inner.y = y + height;
		mark_inner.is_boundary = true;
		(*inner_mark)[height - 1].push_back(mark_inner);
	}
	//MARK mark_right;
	mark_right.x = x + width;
	mark_right.y = y + height;
	mark_right.is_boundary = true;
	(*inner_mark)[height - 1].push_back(mark_right);

}

void ScanningLine::RestScannigLine()
{
	for (int i = 0; i < inner_mark->size(); i++) {
		inner_mark->at(i).clear();
	}
	inner_mark->clear();
	y_max = y_min = 0;
}

void ScanningLine::UpdateAndResortAet(std::vector<EDGE>& aet)
{
	for (int i = 0; i < aet.size(); i++) {
		UpdateAetEdgeInfo(aet[i]);
	}
	std::sort(aet.begin(), aet.end(), MyEdgeXiComparator);
}

void ScanningLine::RemoveNonActiveEdgeFromAet(std::vector<EDGE>& aet, int y)
{
	int i = 0;
	while (i < aet.size()) {
		if (MyIsEdgeOutOfActive(aet[i], y)) {
			vector<EDGE>::iterator it = aet.begin();
			for (int j = 0; j < i; j++) {
				it++;
			}
			aet.erase(it);
		}
		else {
			i++;
		}
	}
}