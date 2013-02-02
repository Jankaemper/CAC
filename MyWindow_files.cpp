//Methods handling files (open, load, save)

// Copyright (c) 2005  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: svn+ssh://scm.gforge.inria.fr/svn/cgal/trunk/Arrangement_on_surface_2/demo/Arrangement_on_surface_2/MyWindow_files.cpp $
// $Id: MyWindow_files.cpp 57467 2010-07-12 10:02:56Z glisse $
//
//
//
// Author(s)     : Baruch Zukerman <baruchzu@post.tau.ac.il>

#include <CGAL/basic.h>


#include "arrangement_2.h"
#include "forms.h"
#include "qt_layer.h"
#include "demo_tab.h"
#include <string>
using namespace std;


#ifdef CGAL_USE_CORE
#include "Conic_reader.h"
#endif

#include <CGAL/IO/Arr_with_history_iostream.h>

/*! open a segment file and add new tab*/
void MyWindow::fileOpenSegment()
{
  Qt_widget_base_tab *w_demo_p =
    static_cast<Qt_widget_base_tab *> (myBar->currentPage());
  bool flag = (w_demo_p->traits_type == CONIC_TRAITS);
  if( w_demo_p->is_empty() ) // pm is empty
  {
    fileOpen(true,true);
  }
  else
  {
    FileOpenOptionsForm
         *form = new FileOpenOptionsForm(flag);
    if ( form->exec() )
    {
        int id = form->buttonGroup->id(form->buttonGroup->selected());
		switch ( id )
		{
          case 0: // open file in a new tab
            add_conic_tab();
            fileOpen(true,false);
            break;
          case 1: // open file in current tab (delete current Pm)
            fileOpen(true,true);
          break;
          case 2: // merge file into current tab
            fileOpen(true,false);
          break;
		}// switch
	}// if

  }

}// fileOpenSegment

/*! open a segment file and add new tab*/
void MyWindow::fileOpenConic()
{

  Qt_widget_base_tab *w_demo_p =
    static_cast<Qt_widget_base_tab *> (myBar->currentPage());
  bool flag = (w_demo_p->traits_type == CONIC_TRAITS);
  if( w_demo_p->is_empty() ) // pm is empty
  {
    fileOpen(false, true);
  }
  else
  {
    FileOpenOptionsForm
         *form = new FileOpenOptionsForm(flag);
    if ( form->exec() )
    {
        int id = form->buttonGroup->id(form->buttonGroup->selected());
		switch ( id )
		{
          case 0: // open file in a new tab
            add_conic_tab();
            fileOpen(false,false);
            break;
          case 1: // open file in current tab (delete current Pm)
            fileOpen(false, true);
          break;
          case 2: // merge file into current tab
            fileOpen(false,false);
          break;
		}// switch
	}// if
  }

}// fileOpenSegment


/*! open a file */
void MyWindow::fileOpen( bool indicator,bool clear_flag )
{
 printf("enter file open \n");
  QString filename =
    QFileDialog::getOpenFileName(QString::null, 0, this,
                                 "file open", "Demo -- File Open" );
  if ( !filename.isEmpty() )
  {
    load(indicator, filename , clear_flag);
	updateMode( dragMode );
  }
  else
    statusBar()->message( "File Open abandoned", 2000 );
  printf("quit file open \n");

}

/*! open a Pm file */
void MyWindow::fileOpenPm()
{
  QString filename =
    QFileDialog::getOpenFileName(QString::null, 0, this,
                                 "file open", "Demo -- File Open" );
  if ( filename.isEmpty() )
  {
    statusBar()->message( "File Open abandoned", 2000 );
	return;
  }
  std::ifstream inputFile(filename.ascii());
  // Creates an ifstream object named inputFile
  if (! inputFile.is_open()) // Always test file open
  {
    std::cout << "Error opening input file" << std::endl;
    return;
  }

  Qt_widget_base_tab    *w_demo_p1 =
    static_cast<Qt_widget_base_tab *> (myBar->currentPage());

  QCursor old = w_demo_p1->cursor();
  w_demo_p1->setCursor(Qt::WaitCursor);
  w_demo_p1->read_from_file = true;
  switch ( w_demo_p1->traits_type ) {
   case SEGMENT_TRAITS:
    {
      Qt_widget_demo_tab<Segment_tab_traits> *w_demo_p =
       static_cast<Qt_widget_demo_tab<Segment_tab_traits> *>
       (myBar->currentPage());
      inputFile >> (*w_demo_p->m_curves_arr);
     break;
    }
   case POLYLINE_TRAITS: // dosen't work !!
    {
      Qt_widget_demo_tab<Polyline_tab_traits> *w_demo_p =
       static_cast<Qt_widget_demo_tab<Polyline_tab_traits> *>
       (myBar->currentPage());
     inputFile >> (*w_demo_p->m_curves_arr);
	   break;
    }
   case CONIC_TRAITS: // dosen't work !!
    {
      /*Qt_widget_demo_tab<Conic_tab_traits> *w_demo_p =
       static_cast<Qt_widget_demo_tab<Conic_tab_traits> *>
       (myBar->currentPage());
     inputFile >> (*w_demo_p->m_curves_arr);*/
     break;
    }
  }

  inputFile.close();


  w_demo_p1->set_window(w_demo_p1->bbox.xmin() , w_demo_p1->bbox.xmax() ,
                     w_demo_p1->bbox.ymin() , w_demo_p1->bbox.ymax());

  inputFile.close();
  w_demo_p1->setCursor(old);
  updateMode( dragMode );
  something_changed();

  setCaption( QString( "Arrangement -- %1" ).arg( m_filename ) );
  statusBar()->message( QString( "Opened \'%1\'" ).arg( m_filename ), 2000 );

}
/*! open a polyline or conic file
 * \param filename - name of the file
 */
void MyWindow::load( bool indicator,const QString& filename , bool clear_flag )
{
  // bring into correct format
  QString filenameNew;
  string stringName = filename.ascii();
  if (indicator == true && stringName.find("Zerlegt")==-1)
  {
	  filenameNew= QString(prescanInput(filename));
  }
  else
  {
	  filenameNew= QString(filename);
  }
  std::ifstream inputFile(filenameNew.ascii());
  // Creates an ofstream object named inputFile
  if (! inputFile.is_open()) // Always test file open
  {
    std::cout << "Error opening input file" << std::endl;
    return;
  }

  Qt_widget_base_tab    *w_demo =
    static_cast<Qt_widget_base_tab *> (myBar->currentPage());

  QCursor old = w_demo->cursor();
  w_demo->setCursor(Qt::WaitCursor);

  if (indicator==true)
  {
    Qt_widget_demo_tab<Conic_tab_traits>    *w_demo_p =
      static_cast<Qt_widget_demo_tab<Conic_tab_traits> *>
      (myBar->currentPage());
    if (clear_flag)
        w_demo_p->m_curves_arr->clear();

    int count;
    inputFile >> count;
    int i;

    std::list<Arr_conic_2> seg_list;
    for (i = 0; i < count; i++) {
      float f0, f1, f2, f3;
      inputFile >> f0 >> f1 >> f2 >> f3;

      Rational x0 (static_cast<int> (10000*f0+0.5), 10000);
      Rational y0 (static_cast<int> (10000*f1+0.5), 10000);
      Rational x1 (static_cast<int> (10000*f2+0.5), 10000);
      Rational y1 (static_cast<int> (10000*f3+0.5), 10000);

      Rat_point_2 p1(x0,y0);
      Rat_point_2 p2(x1,y1);

      Rat_segment_2 seg(p1, p2);
	  Arr_conic_2 curve(seg);


      CGAL::Bbox_2 curve_bbox = curve.bbox();
      if (i == 0)
        w_demo->bbox = curve_bbox;
      else
        w_demo->bbox = w_demo->bbox + curve_bbox;

      seg_list.push_back(curve);
    }

    CGAL::insert(*(w_demo_p->m_curves_arr), seg_list.begin(), seg_list.end());
  }
  else
  {
    Qt_widget_demo_tab<Conic_tab_traits>    *w_demo_p =
      static_cast<Qt_widget_demo_tab<Conic_tab_traits> *>
      (myBar->currentPage());
    if (clear_flag)
        w_demo_p->m_curves_arr->clear();
    typedef Conic_tab_traits::Traits    Conic_traits;
    Conic_reader<Conic_traits>  reader;
    std::list<Arr_conic_2>               curve_list;
    reader.read_data(filename, std::back_inserter(curve_list), w_demo->bbox);
    CGAL::insert (*(w_demo_p->m_curves_arr), curve_list.begin(),
                  curve_list.end());
  }

  w_demo->set_window(w_demo->bbox.xmin() , w_demo->bbox.xmax() ,
                     w_demo->bbox.ymin() , w_demo->bbox.ymax());

  inputFile.close();
  w_demo->setCursor(old);

  something_changed();

}


QString MyWindow::prescanInput( const QString& in_filename )
{

  std::ifstream inputFile(in_filename.ascii());
  if (! inputFile.is_open())
  {
    std::cout << "Error opening input file" << std::endl;
    return in_filename;
  }

  QString out_filename (in_filename.ascii());
  QString stringAdd("erlegt");
  string out_filename_string =out_filename.ascii();
  int posInsert = out_filename_string.find("/Z_");
  out_filename.insert(posInsert+2,stringAdd);
  std::ofstream outputFile(out_filename.ascii());

  string x1 = "";
  string y1 = "";
  string x2 = "";
  string y2 = "";
  string xm = "";
  string ym = "";
  string line;
  int linesCounter = 0;
  while (inputFile.good())
  {
	  getline(inputFile,line);
	  if (line[0]=='<')
	  {
		  x1 = "";
		  x2 = "";
		  y1 = "";
		  y2 = "";
	  }
	  else if (line[0]=='h')
	  {
		  outputFile << xm<< ym << x2 << y2 << "\n";
		  linesCounter++;
	  }
	  else if (x1=="")
	  {
		  int firstSpace = line.find(" ");
		  x1 = line.substr(0,firstSpace+1);
		  xm = line.substr(0,firstSpace+1);
		  int secSpace = line.find(" ",firstSpace+1);
		  y1 = line.substr(firstSpace+1,secSpace-firstSpace);
		  ym = line.substr(firstSpace+1,secSpace-firstSpace);
	  }
	  else if (x2 =="")
	  {
		  int firstSpace = line.find(" ");
		  x2 = line.substr(0,firstSpace+1);
		  int secSpace = line.find(" ",firstSpace+1);
		  y2 = line.substr(firstSpace+1,secSpace-firstSpace);
		  outputFile << x1 << y1 << x2 << y2 << "\n";
		  linesCounter++;
	  }
	  else
	  {
		  x1 = x2;
		  y1 = y2;
		  int firstSpace = line.find(" ");
		  x2 = line.substr(0,firstSpace+1);
		  int secSpace = line.find(" ",firstSpace+1);
		  y2 = line.substr(firstSpace+1,secSpace-firstSpace);
		  outputFile << x1 << y1 << x2 << y2 << "\n";
		  linesCounter++;

	  }
  }


  //append the number of lines to the beginning of output file


  inputFile.close();
  outputFile.close();
  return out_filename;
}


#ifdef CGAL_USE_CORE
/*! read a conic curve
 * \param is - input file stream
 * \param cv - will hold the reading curve
 */
void MyWindow::ReadCurve(std::ifstream & is, Arr_conic_2 & cv)
{
  // Read a line from the input file.
  char one_line[128];

  skip_comments (is, one_line);
  std::istringstream str_line (one_line);

  // Read the arc type and act accordingly.
  char     type;

  str_line >> type;

  if (type == 's' || type == 'S')
  {
    // Construct a line segment. The line should have the format:
    //   s <x1> <y1> <x2> <y2>
    // where (x1, y1), (x2, y2) are the endpoints of a segment.
    Rational    x1, y1, x2, y2;

    str_line >> x1 >> y1 >> x2 >> y2;

    Rat_point_2   p1(x1, y1), p2(x2, y2);
    Rat_segment_2 seg (p1, p2);

    cv = Arr_conic_2 (seg);
  }
  else if (type == 'c' || type == 'C')
  {
    // Construct a full circle. The line should have the format:
    //   c <x0> <y0> <R_sq>
    // where (x0, y0) is the center of the circle and R_sq is its squared
    // radius.
    Rational    x0, y0, R_sq;
    str_line >> x0 >> y0 >> R_sq;

    Rat_point_2   p0(x0, y0);
    Rat_circle_2  circ(p0, R_sq);

    cv = Arr_conic_2 (circ);
  }
  else if (type == 't' || type == 'T')
  {
    // Construct a circular arc. The line should have the format:
    //   t <x1> <y1> <x2> <y2> <x3> <y3>
    // where (x1, y1), (x2, y2) and (x3, y3) define the arc.
    Rational    x1, y1, x2, y2, x3, y3;

    str_line >> x1 >> y1 >> x2 >> y2 >> x3 >> y3;

    Rat_point_2   p1(x1, y1), p2(x2, y2), p3(x3, y3);

    cv = Arr_conic_2 (p1, p2, p3);
  }
  else if (type == 'f' || type == 'F')
  {
    // Construct a full conic curve. The line should have the format:
    //   c <r> <s> <t> <u> <v> <w>
    // where r, s, t, u, v, w define the conic equation.
    Rational    r, s, t, u, v, w;

    str_line >> r >> s >> t >> u >> v >> w;

    cv = Arr_conic_2 (r, s, t, u, v, w);
  }
  else if (type == 'a' || type == 'A')
  {
    // Construct a conic arc. The line should have the format:
    //   c <r> <s> <t> <u> <v> <w> <orient> <x1> <y1> <x2> <y2>
    // where r, s, t, u, v, w define the conic equation, while (x1, y1)
    // and (x2, y2) are the arc's endpoints.
    Rational    r, s, t, u, v, w;

    str_line >> r >> s >> t >> u >> v >> w;

    // Read the orientation.
    int               i_orient;
    CGAL::Orientation orient;

    str_line >> i_orient;
    if (i_orient > 0)
      orient = CGAL::COUNTERCLOCKWISE;
    else if (i_orient < 0)
      orient = CGAL::CLOCKWISE;
    else
      orient = CGAL::COLLINEAR;

    // Read the end points of the arc and create it.
    // Notice we read the coordinates as strings, then we convert them to
    // the Algebraic type, as we do not want to initialize Algebraic from
    // a double.
    char    num[50];
    Algebraic    x1, y1, x2, y2;

    str_line >> num;
    x1 = Algebraic(num);
    str_line >> num;
    y1 = Algebraic(num);

    str_line >> num;
    x2 = Algebraic(num);
    str_line >> num;
    y2 = Algebraic(num);

    Arr_conic_point_2 ps (x1, y1);
    Arr_conic_point_2 pt (x2, y2);

    cv = Arr_conic_2 (r, s, t, u, v, w, orient, ps ,pt);
  }
  else if (type == 'l' || type == 'L')
  {
    // Construct a conic arc. The line should have the format:
    //   c <r> <s> <t> <u> <v> <w> <a> <b> <c>
    // where r, s, t, u, v, w define the conic equation and a, b, c define
    // a line that intersects it.
    Rational    r, s, t, u, v, w;
    Rational    a, b, c;

    str_line >> r >> s >> t >> u >> v >> w >> a >> b >> c;

    Rat_line_2    line (a, b, c);

  }
  else if (type == 'q' || type == 'Q')
  {
    // Construct a circular arc. The line should have the format:
    //   t <x1> <y1> <x2> <y2> <x3> <y3> <x4> <y4> <x5> <y5>
    // where (x1, y1), (x2, y2), (x3, y3), (x4, y4) and (x5, y5) define the
    // arc.
    Rational    x1, y1, x2, y2, x3, y3, x4, y4, x5, y5;

    str_line >> x1 >> y1 >> x2 >> y2 >> x3 >> y3 >> x4 >> y4 >> x5 >> y5;

    Rat_point_2   p1(x1, y1), p2(x2, y2), p3(x3, y3), p4(x4, y4), p5(x5, y5);

    cv = Arr_conic_2 (p1, p2, p3, p4, p5);
  }
  else
  {
    std::cerr << "Illegal conic type specification: " << type << "."
	      << std::endl;
  }

  return;
}
#endif

/*! calls from readCurve to skip comments in the input file */
void MyWindow::skip_comments( std::ifstream& is, char* one_line )
{
  while( !is.eof() )
  {
    is.getline( one_line, 128 );
    if( one_line[0] != '#' )
    {
      break;
    }
  }
}

/*! save file in a different name */
void MyWindow::fileSaveAs()
{
  QString filename =
    QFileDialog::getSaveFileName(QString::null, "Arrangement (*.arr)", this,
                                 "file save as", "Arrangement -- File Save As");
  if ( !filename.isEmpty() )
  {
    int answer = 0;
    if ( QFile::exists( filename ) )
      answer = QMessageBox::warning(
                                    this, "Overwrite File",
                                    QString( "Overwrite\n\'%1\'?" ).
                                    arg( filename ),
                                    "&Yes", "&No", QString::null, 1, 1 );
    if ( answer == 0 )
    {
      m_filename = filename;
      //updateRecentFiles( filename );
      fileSave();
      return;
    }
  }
  statusBar()->message( "Saving abandoned", 2000 );
}

/*! save file */
void MyWindow::fileSave()
{
  if ( m_filename.isEmpty() ) {
    fileSaveAs();
    return;
  }

  std::ofstream outFile(m_filename.ascii());
  // Creates an ofstream object named outFile
  if (! outFile.is_open()) // Always test file open
  {
    std::cout << "Error opening input file" << std::endl;
    return;
  }

  Qt_widget_base_tab    *w_demo_p1 =
    static_cast<Qt_widget_base_tab *> (myBar->currentPage());

  switch ( w_demo_p1->traits_type ) {
   case SEGMENT_TRAITS:
    {
     Qt_widget_demo_tab<Segment_tab_traits> *w_demo_p =
       static_cast<Qt_widget_demo_tab<Segment_tab_traits> *>
       (myBar->currentPage());
     outFile << (*w_demo_p->m_curves_arr);
     break;
    }
   case POLYLINE_TRAITS:
    {
      Qt_widget_demo_tab<Polyline_tab_traits> *w_demo_p =
       static_cast<Qt_widget_demo_tab<Polyline_tab_traits> *>
       (myBar->currentPage());
      outFile << (*w_demo_p->m_curves_arr);
     break;
    }
   case CONIC_TRAITS:
    {
      /*Qt_widget_demo_tab<Conic_tab_traits> *w_demo_p =
       static_cast<Qt_widget_demo_tab<Conic_tab_traits> *>
       (myBar->currentPage());
      outFile << (*w_demo_p->m_curves_arr);*/
     break;
    }
  }

  outFile.close();

  setCaption( QString( "Arrangement -- %1" ).arg( m_filename ) );
  statusBar()->message( QString( "Saved \'%1\'" ).arg( m_filename ), 2000 );
}

/*! save planar map to post script */
void MyWindow::fileSave_ps()
{
#if 0

  Qt_widget_base_tab    *w_demo_p1 =
    static_cast<Qt_widget_base_tab *> (myBar->currentPage());

  switch ( w_demo_p1->traits_type ) {
   case SEGMENT_TRAITS:
    {
     Qt_widget_demo_tab<Segment_tab_traits> *w_demo_p =
       static_cast<Qt_widget_demo_tab<Segment_tab_traits> *>
       (myBar->currentPage());
     // Print to Postscript file:
     CGAL::Postscript_file_stream  LPF(m_width, m_height ,"pm.ps");
     LPF.init(-3,3,-3);
     LPF.set_line_width(1);
     //LPF << w_demo_p->m_curves_arr;
     break;
    }
   case POLYLINE_TRAITS:
    {
     //Qt_widget_demo_tab<Polyline_tab_traits> *w_demo_p =
     //static_cast<Qt_widget_demo_tab<Polyline_tab_traits> *>
     //  (myBar->currentPage());
     //outFile << w_demo_p->m_curves_arr;
     break;
    }
   case CONIC_TRAITS:
    {
     //Qt_widget_demo_tab<Conic_tab_traits> *w_demo_p =
     //static_cast<Qt_widget_demo_tab<Conic_tab_traits> *>
      // (myBar->currentPage());
     //outFile << w_demo_p->m_curves_arr;
     break;
    }
  }
#endif
}

