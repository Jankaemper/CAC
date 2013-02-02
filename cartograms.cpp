
#include <CGAL/basic.h>
#include "arrangement_2.h"
#include "forms.h"
#include "qt_layer.h"
#include "demo_tab.h"
#include <CGAL/Polygon_2.h>
#include <CGAL/IO/Qt_widget_Polygon_2.h>
#include "cgal_types.h"
#include <CGAL/Arrangement_2.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arr_extended_dcel.h>
#include <CGAL/Sweep_line_2_algorithms.h>
#include <list>
#include <CGAL/Arrangement_with_history_2.h>
#include <math.h>
#include "print.h"
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <string>
#include <map>
#include <sys/time.h>
using namespace std;


//straight skeleton includes
#include<boost/shared_ptr.hpp>
#include<CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include<CGAL/create_straight_skeleton_2.h>
#include<CGAL/create_straight_skeleton_from_polygon_with_holes_2.h>
#include<CGAL/Polygon_with_holes_2.h>



//dual graph in boost includes and flow algo includes


#include <climits>
#include <boost/graph/graph_utility.hpp>
#include <dual_history_bug.h>
#include <utility>
#include <vector>
#include <algorithm> // for std::min and std::max
#include <functional>
#include <Arr_face_index_map.h>
#include <CGAL/graph_traits_Dual_Arrangement_2.h>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/config.hpp>
#include <boost/bind.hpp>
#include <boost/graph/vector_as_graph.hpp>
#include <boost/graph/strong_components.hpp>
#include <boost/graph/topological_sort.hpp>
#include <boost/graph/named_function_params.hpp>
#include <boost/graph/push_relabel_max_flow.hpp>
#include <boost/graph/boykov_kolmogorov_max_flow.hpp>
#include <boost/graph/graph_concepts.hpp>
#include <boost/graph/read_dimacs.hpp>
#include <boost/graph/adjacency_list.hpp>



#define PI 3.14159265

/*
 creates the dual graph of the arrangement/subdivision, adds source and sink and reverse edges (required by boost
 algorithms) and writes capacities
*/

dualGraph MyWindow::extractDualFromArrangement(primalMap *indPrimalFaces, dualMap *indDualVertices)
{
	//dual graph can is extracted manually from arrangement (iterating over faces and outer boundaries) and then written into a adjacency list
	//weights, names, capacities are attached with property maps
	//the flow is implicitly obtained from a property map also (residual capacity)

	Qt_widget_demo_tab<Conic_tab_traits> *panel =
		static_cast<Qt_widget_demo_tab<Conic_tab_traits> *> (myBar->currentPage());

	Conic_arr arr (*panel->m_curves_arr);
	Conic_face_iterator fit;
	Conic_face_handle fhan;
	int i = 0;

	primalMap indexPrimalFaces = *indPrimalFaces;

	dual_edge dummy;

	int number_of_faces = arr.number_of_faces();
	vector<vector<dual_edge> > adjMatrix;

	for ( i = 0; i < number_of_faces; i++ )
	{
		adjMatrix.push_back ( vector<dual_edge>() );
		for ( int j = 0; j < number_of_faces; j++ )
		{
			adjMatrix[i].push_back ( dummy );
		}
	}

	dualGraph dual;
	dual_vertex dualVertex;
	dual_vertex outer2;
	dual_edge dualEdge, dualEdge2;
	primal_edge primalEdge;

	//get property maps for vertex and edge properties
	property_map<dualGraph, edge_capacity_t>::type capacity = get(edge_capacity_t(), dual);
	property_map<dualGraph, vertex_name_t>::type name = get(vertex_name_t(), dual);
	property_map<dualGraph, vertex_weight_t>::type weight = get(vertex_weight_t(), dual);
	property_map<dualGraph, vertex_initArea_t>::type initArea = get(vertex_initArea_t(), dual);
	property_map<dualGraph, vertex_conti_t>::type conti = get(vertex_conti_t(), dual);

	property_map<dualGraph, edge_residual_capacity_t>::type res_capacity = get(edge_residual_capacity_t(), dual);
	property_map<dualGraph, edge_reverse_t>::type reverse = get(edge_reverse_t(), dual);
	property_map<dualGraph, edge_primalEdges_t>::type primals = get(edge_primalEdges_t(), dual);
	property_map<dualGraph, edge_primalEdgesCap_t>::type primalCaps = get(edge_primalEdgesCap_t(), dual);
	property_map<dualGraph, edge_primalEdgesSources_t>::type primalSources = get(edge_primalEdgesSources_t(), dual);
	property_map<dualGraph, edge_primalEdgesTargets_t>::type primalTargets = get(edge_primalEdgesTargets_t(), dual);
	property_map<dualGraph, edge_primalEdgesFlow_t>::type primalFlow = get(edge_primalEdgesFlow_t(), dual);


	// create mapping of indices to dual vertices with the correct indices (given in primal arrangement)
	// add all the vertices to the dual graph (iterate over arrangement faces)
	for (fit = panel->m_curves_arr->faces_begin(); fit != panel->m_curves_arr->faces_end(); ++fit)
	{
		dualVertex = boost::add_vertex(dual);
		(*indDualVertices)[indexPrimalFaces[fit]] = dualVertex;

		//set the properties of the dual vertex (see cgal_types for definition)
		put(name, dualVertex, fit->getCountry());
		put(weight, dualVertex, fit->getWeight());
		put(initArea,dualVertex,fit->getInitialArea());
		put(conti,dualVertex,fit->getConti());


		//printf("added dual vertex for %s with index %i | weight = %f | area = %f \n", fit->getCountry(), indexPrimalFaces[fit],get(weight,dualVertex),get(initArea,dualVertex));
	}

	//search for existing original outer face vertex
	dual_vertex outer1;
	boost::graph_traits<dualGraph>::vertex_iterator vi, vi_end, next;
	boost::tie(vi, vi_end) = boost::vertices(dual);
	for (next = vi; vi != vi_end; vi = next)
	{
		++next;
		if (get(name,*vi) == "OUTERFACE")
		{
			outer1 = *vi;
		}
	}

	//add second outer face vertex
	outer2 = boost::add_vertex(dual);
	put(name, outer2,"OUTERFACE2");
	put(weight, outer2,0);
	put(initArea,outer2,0);
	put(conti,dualVertex,-1);




	// add all the edges to the dual graph (iterate over arrangement faces and then iterate over outer boundary)
	for (fit = panel->m_curves_arr->faces_begin(); fit != panel->m_curves_arr->faces_end(); ++fit)
	{
		if (!fit->is_unbounded())
		{
			Conic_ccb_halfedge_circulator cc=fit->outer_ccb();
			//printf("\n \n Edges from %s \n",fit->getCountry());
			do {
				  //discard edges that lie completely in the face (antennas)
				  if (panel->antenna(cc))
					continue;
				  Conic_halfedge_handle hhan = cc;

				  //get source and target index and compare (if target index is smaller, then this edge pair has already been added in a previous iteration)
				  //so in that case we skip the rest

				  if (indexPrimalFaces[fit] < indexPrimalFaces[hhan->twin()->face()] || hhan->twin()->face()->getCountry() == "OUTERFACE")
				  {
					  // we dont allow parallel edges in the dual, so look up in the adjacency matrix
					  // if a new pair has to be created or has already been created in previous iteration
					  if (adjMatrix[indexPrimalFaces[fit]][indexPrimalFaces[hhan->twin()->face()]] == dummy)
					  {
						  if (hhan->twin()->face()->getCountry() == "OUTERFACE")
						  {
							  //case of outer face with special treatment
							  //we have two vertices representing the outer face: one that is connected with growing faces and one that is connected with shrinking faces
							  if (get(weight,(*indDualVertices)[indexPrimalFaces[fit]]) >= get(initArea,(*indDualVertices)[indexPrimalFaces[fit]]))
							  {
								  //face wants to grow: connect with outer1 (the original outer face vertex)
								  dualEdge = (boost::add_edge((*indDualVertices)[indexPrimalFaces[fit]],(*indDualVertices)[indexPrimalFaces[hhan->twin()->face()]],dual)).first;
								  put(capacity, dualEdge, calcBendingCapacityExact(hhan->twin()->face(),hhan->twin(),indPrimalFaces));
								  dualEdge2 = (boost::add_edge((*indDualVertices)[indexPrimalFaces[hhan->twin()->face()]],(*indDualVertices)[indexPrimalFaces[fit]],dual)).first;
								  put(capacity, dualEdge2, calcBendingCapacityExact(hhan->face(),hhan,indPrimalFaces));
								  put(reverse,dualEdge,dualEdge2);
								  put(reverse,dualEdge2,dualEdge);
								  adjMatrix[indexPrimalFaces[fit]][indexPrimalFaces[hhan->twin()->face()]] = dualEdge;
								  get(primals,dualEdge).push_back(hhan);
								  get(primalCaps,dualEdge).push_back(calcBendingCapacityExact(hhan->twin()->face(),hhan->twin(),indPrimalFaces));
								  get(primals,dualEdge2).push_back(hhan->twin());
								  get(primalCaps,dualEdge2).push_back(calcBendingCapacityExact(hhan->face(),hhan,indPrimalFaces));
								  //store the vertices seperatly and redondantly (needed for storage of svg graphic after bending)
								  get(primalSources,dualEdge).push_back(Arr_conic_point_2(hhan->source()->point()));
								  get(primalTargets,dualEdge).push_back(Arr_conic_point_2(hhan->target()->point()));
								  get(primalSources,dualEdge2).push_back(Arr_conic_point_2(hhan->target()->point()));
								  get(primalTargets,dualEdge2).push_back(Arr_conic_point_2(hhan->source()->point()));


							  }
							  else
							  {
								  //face wants to shrink: connect with outer2 (the auxiliary outer face vertex)
								  dualEdge = (boost::add_edge((*indDualVertices)[indexPrimalFaces[fit]],outer2,dual)).first;
								  put(capacity, dualEdge, calcBendingCapacityExact(hhan->twin()->face(),hhan->twin(),indPrimalFaces));
								  dualEdge2 = (boost::add_edge(outer2,(*indDualVertices)[indexPrimalFaces[fit]],dual)).first;
								  put(capacity, dualEdge2, calcBendingCapacityExact(hhan->face(),hhan,indPrimalFaces));
								  put(reverse,dualEdge,dualEdge2);
								  put(reverse,dualEdge2,dualEdge);
								  adjMatrix[indexPrimalFaces[fit]][indexPrimalFaces[hhan->twin()->face()]] = dualEdge;
								  get(primals,dualEdge).push_back(hhan);
								  get(primalCaps,dualEdge).push_back(calcBendingCapacityExact(hhan->twin()->face(),hhan->twin(),indPrimalFaces));
								  get(primals,dualEdge2).push_back(hhan->twin());
								  get(primalCaps,dualEdge2).push_back(calcBendingCapacityExact(hhan->face(),hhan,indPrimalFaces));
								  //store the vertices seperatly and redondantly (needed for storage of svg graphic after bending)
								  get(primalSources,dualEdge).push_back(Arr_conic_point_2(hhan->source()->point()));
								  get(primalTargets,dualEdge).push_back(Arr_conic_point_2(hhan->target()->point()));
								  get(primalSources,dualEdge2).push_back(Arr_conic_point_2(hhan->target()->point()));
								  get(primalTargets,dualEdge2).push_back(Arr_conic_point_2(hhan->source()->point()));
							  }
						  }
						  else
						  {
							  //both faces are bounded
							  dualEdge = (boost::add_edge((*indDualVertices)[indexPrimalFaces[fit]],(*indDualVertices)[indexPrimalFaces[hhan->twin()->face()]],dual)).first;
							  //store capacity and cost as bundled edge property
							  put(capacity, dualEdge, calcBendingCapacityExact(hhan->twin()->face(),hhan->twin(),indPrimalFaces));
							  //do the same for reverse edge
							  dualEdge2 = (boost::add_edge((*indDualVertices)[indexPrimalFaces[hhan->twin()->face()]],(*indDualVertices)[indexPrimalFaces[fit]],dual)).first;
							  put(capacity, dualEdge2, calcBendingCapacityExact(hhan->face(),hhan,indPrimalFaces));

							  //set the reverse edges accordingly
							  put(reverse,dualEdge,dualEdge2);
							  put(reverse,dualEdge2,dualEdge);
							  adjMatrix[indexPrimalFaces[fit]][indexPrimalFaces[hhan->twin()->face()]] = dualEdge;
							  //printf("added edge to %s with capacity %f \n",hhan->twin()->face()->getCountry(), get(capacity,dualEdge) );
							  //printf("added edge from %s with capacity %f \n",hhan->twin()->face()->getCountry(), get(capacity,dualEdge2) );

							  //store primal arrangement edges that are associated with current dual edge and also store their single capacities
							  get(primals,dualEdge).push_back(hhan);
							  get(primalCaps,dualEdge).push_back(calcBendingCapacityExact(hhan->twin()->face(),hhan->twin(),indPrimalFaces));
							  get(primals,dualEdge2).push_back(hhan->twin());
							  get(primalCaps,dualEdge2).push_back(calcBendingCapacityExact(hhan->face(),hhan,indPrimalFaces));

							  //store the vertices seperatly and redondantly (needed for storage of svg graphic after bending)
							  get(primalSources,dualEdge).push_back(Arr_conic_point_2(hhan->source()->point()));
							  get(primalTargets,dualEdge).push_back(Arr_conic_point_2(hhan->target()->point()));
							  get(primalSources,dualEdge2).push_back(Arr_conic_point_2(hhan->target()->point()));
							  get(primalTargets,dualEdge2).push_back(Arr_conic_point_2(hhan->source()->point()));

						  }

					  }
					  //if the pair already exists we increase the capacity of the edge
					  else
					  {
						  //cases where both faces are bounded or one is unbounded both
						  dualEdge  = adjMatrix[indexPrimalFaces[fit]][indexPrimalFaces[hhan->twin()->face()]];
						  put(capacity, dualEdge, get(capacity,dualEdge) + calcBendingCapacityExact(hhan->twin()->face(),hhan->twin(),indPrimalFaces));
						  dualEdge2  = get(reverse,dualEdge);
						  put(capacity, dualEdge2, get(capacity,dualEdge2) + calcBendingCapacityExact(hhan->face(),hhan,indPrimalFaces));
						  //store primal arrangement edges that are associated with current dual edge and also store their single capacities
						  get(primals,dualEdge).push_back(hhan);
						  get(primalCaps,dualEdge).push_back(calcBendingCapacityExact(hhan->twin()->face(),hhan->twin(),indPrimalFaces));
						  get(primals,dualEdge2).push_back(hhan->twin());
						  get(primalCaps,dualEdge2).push_back(calcBendingCapacityExact(hhan->face(),hhan,indPrimalFaces));

						  //store the vertices seperatly and redondantly (needed for storage of svg graphic after bending)
						  get(primalSources,dualEdge).push_back(Arr_conic_point_2(hhan->source()->point()));
						  get(primalTargets,dualEdge).push_back(Arr_conic_point_2(hhan->target()->point()));
						  get(primalSources,dualEdge2).push_back(Arr_conic_point_2(hhan->target()->point()));
						  get(primalTargets,dualEdge2).push_back(Arr_conic_point_2(hhan->source()->point()));

					  }
				  }
			} while (++cc != fit->outer_ccb());
		}
	}

	//now add source and target and connect them to the vertices in the graph that want to grow (source connection) or shrink (target connection)
	dual_vertex source = boost::add_vertex(dual);
	put(name,source,"SOURCE");
	put(weight, source, 0);
	put(initArea,source,0);
	put(conti,dualVertex,-1);

	dual_vertex target = boost::add_vertex(dual);
	put(name,target,"TARGET");
	put(weight, target, 0);
	put(initArea,target,0);
	put(conti,dualVertex,-1);



	//add edge from source to auxiliary outer face (virtual flow)
	dualEdge = (boost::add_edge(outer1, target, dual)).first;
	put(capacity, dualEdge,10000000);
	//we have to add the opposite edge as well (required by push relabel)
	dualEdge2 = (boost::add_edge(target, outer1, dual)).first;
	put(capacity, dualEdge2, 0);
	//set the reverse edges accordingly
	put(reverse,dualEdge,dualEdge2);
	put(reverse,dualEdge2,dualEdge);

	//add edge from original outer face to target (virtual flow)
	dualEdge = (boost::add_edge(source, outer2, dual)).first;
	put(capacity, dualEdge,10000000);
	//we have to add the opposite edge as well (required by push relabel)
	dualEdge2 = (boost::add_edge(outer2,source, dual)).first;
	put(capacity, dualEdge2, 0);
	//set the reverse edges accordingly
	put(reverse,dualEdge,dualEdge2);
	put(reverse,dualEdge2,dualEdge);


	//now set edges from/to source/target to/from the regular vertices
	for (fit = panel->m_curves_arr->faces_begin(); fit != panel->m_curves_arr->faces_end(); ++fit)
	{
		if (!fit->is_unbounded())
		{
			dualVertex = (*indDualVertices)[indexPrimalFaces[fit]];
			if ( get(weight,dualVertex) - get(initArea,dualVertex) >= 0)
			{
				//add edge between source vertex and positive weight increase vertices
				dualEdge = (boost::add_edge(source, dualVertex, dual)).first;
				put(capacity, dualEdge, get(weight,dualVertex) - get(initArea,dualVertex));
				//we have to add the opposite edge as well (required by push relabel)
				dualEdge2 = (boost::add_edge(dualVertex, source, dual)).first;
				put(capacity, dualEdge2, 0);
				//set the reverse edges accordingly
				put(reverse,dualEdge,dualEdge2);
				put(reverse,dualEdge2,dualEdge);

			}
			else
			{
				//add edge between negative weight increase vertices and target
				dualEdge = (boost::add_edge(dualVertex, target, dual)).first;
				put(capacity, dualEdge, get(initArea,dualVertex) - get(weight,dualVertex));
				//we have to add the opposite edge as well (required by push relabel)
				dualEdge2 = (boost::add_edge(target, dualVertex, dual)).first;
				put(capacity, dualEdge2, 0);
				//set the reverse edges accordingly
				put(reverse,dualEdge,dualEdge2);
				put(reverse,dualEdge2,dualEdge);
			}

		}
	}


	return dual;

}

//computes a polygon that bounds the hole subdivision (rectangle)
My_polygon_K MyWindow::getBoundingPolygon()
{
	Qt_widget_demo_tab<Conic_tab_traits> *panel = static_cast<Qt_widget_demo_tab<Conic_tab_traits> *> (myBar->currentPage());
	My_polygon_K pgn;
	Conic_halfedge_handle hit;
	double minX = 100000;
	double minY = 100000;
	double maxX = -100000;
	double maxY = -100000;

	for (  hit = panel->m_curves_arr->halfedges_begin();  hit != panel->m_curves_arr->halfedges_end(); hit++ )
	{
		minX = min(minX,CGAL::to_double(hit->source()->point().x()));
		minX = min(minX,CGAL::to_double(hit->target()->point().x()));

		maxX = max(maxX,CGAL::to_double(hit->source()->point().x()));
		maxX = max(maxX,CGAL::to_double(hit->target()->point().x()));

		minY = min(minY,CGAL::to_double(hit->source()->point().y()));
		minY = min(minY,CGAL::to_double(hit->target()->point().y()));

		maxY = min(maxY,CGAL::to_double(hit->source()->point().y()));
		maxY = max(maxY,CGAL::to_double(hit->target()->point().y()));

	}
	pgn.push_back(K_point(maxX + 100000, maxY + 100000));
	pgn.push_back(K_point(minX - 100000, maxY + 100000));
	pgn.push_back(K_point(minX - 100000, minY - 100000));
	pgn.push_back(K_point(maxX + 100000, minY - 100000));

	return pgn;
}


//for a given continent number this function returns the polygon of the outer polygon of the corresponding continent
My_polygon_K MyWindow::getContinentPolygon(int conti)
{
	Qt_widget_demo_tab<Conic_tab_traits> *panel = static_cast<Qt_widget_demo_tab<Conic_tab_traits> *> (myBar->currentPage());
	My_polygon_K pgn;
	Conic_halfedge_handle hit;
	Conic_vertex_handle vhan;
	int continentNumber;

	//find an edge bordering the outer face (to its left) and continent conti (to its right)
	//store this edge's target vertex in vhan
	for (  int i=outerFaceEdges.size()-1;  i>= 0; i--)
	{
		hit = outerFaceEdges[i];
		if (hit->twin()->face()->getConti() == conti)
		{
			vhan = hit->target();
			i = -1;
		}
	}


	Halfedge_around_vertex_circulator hvit;
	Conic_halfedge_handle hhan;
	pgn.push_back(K_point(CGAL::to_double(hit->source()->point().x()),CGAL::to_double(hit->source()->point().y())));
	pgn.push_back(K_point(CGAL::to_double(hit->target()->point().x()),CGAL::to_double(hit->target()->point().y())));


	//iterate over continent's outer vertices until one entire tour has been performed
	while (vhan != hit->source())
	{
		//iterate over halfedges that have the same target vertex
		//among them there must be one which borders the outer face: take the source of that one as next polygon point
		hvit = vhan->incident_halfedges();
		do
		{
			if (hvit->twin()->face()->is_unbounded())
			{
				hhan = hvit;
			}

		} while (++hvit != vhan->incident_halfedges());
		pgn.push_back(K_point(CGAL::to_double(hhan->source()->point().x()),CGAL::to_double(hhan->source()->point().y())));

		vhan = hhan->source();
	}

	return pgn;
}


//assigns numbers to the different continents and writes them to the faces
void MyWindow::clusterContinents()
{
	Qt_widget_demo_tab<Conic_tab_traits> *panel = static_cast<Qt_widget_demo_tab<Conic_tab_traits> *> (myBar->currentPage());
	Conic_face_handle fhan;
	Conic_face_handle fhan2;

	Conic_ccb_halfedge_circulator cc;
	Conic_halfedge_handle hhan;

	int currentConti = 0;
	int thisConti = 0;

	//at the beginning no continents are set (all number are -1)
	for (  fhan = panel->m_curves_arr->faces_begin();  fhan != panel->m_curves_arr->faces_end(); fhan++ )
	{
		fhan->setConti(-1);
	}

	//now for every face we do a dfs which looks if a number has already been assigned to that continent, if not we assign a new number
	for (  fhan = panel->m_curves_arr->faces_begin();  fhan != panel->m_curves_arr->faces_end(); fhan++ )
	{
		if (!fhan->is_unbounded())
		{
			thisConti = searchConti(fhan);
			if (thisConti == -1)
			{
				fhan -> setConti(currentConti);
				currentConti++;
			}
			else
			{
				fhan -> setConti(thisConti);
			}
			for (fhan2 = panel->m_curves_arr->faces_begin();  fhan2 != panel->m_curves_arr->faces_end(); fhan2++ )
			{
				fhan2->set_visited(false);
			}
		}
	}
	numberOfContinents = currentConti;


}

//recursively goes over all faces of a continent and looks if continent number has already been writen to one of them
int MyWindow::searchConti(Conic_face_handle fhan)
{
	Qt_widget_demo_tab<Conic_tab_traits> *panel = static_cast<Qt_widget_demo_tab<Conic_tab_traits> *> (myBar->currentPage());
	int currentConti;
	fhan->set_visited(true);

	if (fhan->getConti() >= 0)
	{
		return fhan->getConti();
	}

    Conic_ccb_halfedge_circulator cc=fhan->outer_ccb();
	do {
	  if (panel->antenna(cc))
		continue;

	  Conic_halfedge_handle hhan = cc;
	  if (!(hhan->twin()->face()->is_unbounded()) && !(hhan->twin()->face()->visited()))
	  {
		  currentConti = searchConti(hhan->twin()->face());
		  if (currentConti >= 0)
		  {
			  return currentConti;
		  }
	  }

	} while (++cc != fhan->outer_ccb());
	return -1;
}


// helper function (for testing)
void MyWindow::cartogram_start()
{
	printf("enter cartogram start \n");

	Qt_widget_demo_tab<Conic_tab_traits> *panel = static_cast<Qt_widget_demo_tab<Conic_tab_traits> *> (myBar->currentPage());
	Conic_face_iterator fit;
	fit = panel->m_curves_arr->faces_begin();
	outerFaceEdges.clear();
	skeletons.clear();

	Conic_halfedge_iterator hit;
	Conic_halfedge_handle hithand;
	//write all edges that border the outer face to the global vector of this class
	for (  hit = panel->m_curves_arr->halfedges_begin();  hit != panel->m_curves_arr->halfedges_end(); hit++ )
	{
		if (hit->face()->is_unbounded())
		{
			hithand = hit;
			outerFaceEdges.push_back(hithand);
		}
	}
	clusterContinents();

    if (!fit->is_unbounded()){
		fit++;
	}


	showSkeleton(fit);

	storeStraight("ouput.ipe");

}


// creates a new tab with a cartogram instance
void MyWindow::cartogram_balance()
{

	printf("enter carto balance \n");

	debug = false;

	Qt_widget_demo_tab<Conic_tab_traits> *panel =
		static_cast<Qt_widget_demo_tab<Conic_tab_traits> *> (myBar->currentPage());
	Conic_arr arr (*panel->m_curves_arr);

	primalMap indexPrimalFaces (*panel->m_curves_arr);
	dualMap indexDualVertices;
	setUpWeights();

	//set up a vector of all edges bordering the outer face in order to avoid multiple requests of this type (no default cgal function for this)
	Conic_halfedge_iterator hit;
	Conic_halfedge_handle hithand;
	outerFaceEdges.clear();
	skeletons.clear();

	//write all edges that border the outer face to the global vector of this class
	for (  hit = panel->m_curves_arr->halfedges_begin();  hit != panel->m_curves_arr->halfedges_end(); hit++ )
	{
		if (hit->face()->is_unbounded())
		{
			hithand = hit;
			outerFaceEdges.push_back(hithand);
		}
	}
	printf(" %i outer face edges \n",(int)outerFaceEdges.size());

	//cluster the continents of the subdivsion (assign number to the faces)
	clusterContinents();
	printf(" %i continents \n",numberOfContinents);

	/*
	Conic_face_handle fit;
	for (  fit = panel->m_curves_arr->faces_begin();  fit != panel->m_curves_arr->faces_end(); fit++ )
	{
		std::cout << fit->getCountry() << " belongs to continent " << fit->getConti() << std::endl;
	}
	 */

	timeval start, end;
	gettimeofday(&start, 0);
	dualGraph dual = extractDualFromArrangement(&indexPrimalFaces,&indexDualVertices);
	gettimeofday(&end, 0);
	  cout << start.tv_sec << ':' << start.tv_usec << endl;
	  cout << end.tv_sec << ':' << end.tv_usec << endl;
	  cout << "duration setting up dual graph: " << end.tv_sec-start.tv_sec << ":" << end.tv_usec-start.tv_usec << endl;


	printf("Size of dual map: %i \n",(int)indexDualVertices.size());

	printDualToFile(dual);

	//get property maps
	property_map<dualGraph, vertex_name_t>::type name = get(vertex_name_t(), dual);
	property_map<dualGraph, edge_capacity_t>::type capacity = get(edge_capacity_t(), dual);
	property_map<dualGraph, edge_residual_capacity_t>::type res_capacity = get(edge_residual_capacity_t(), dual);
	property_map<dualGraph, edge_reverse_t>::type reverse = get(edge_reverse_t(), dual);
	property_map<dualGraph, vertex_weight_t>::type weight = get(vertex_weight_t(), dual);
	property_map<dualGraph, edge_primalEdges_t>::type primals = get(edge_primalEdges_t(), dual);
	property_map<dualGraph, edge_primalEdgesCap_t>::type primalCaps = get(edge_primalEdgesCap_t(), dual);

	double flow;

	dual_vertex source;
	dual_vertex target;
	//get source and target vertex (search by name)
	boost::graph_traits<dualGraph>::vertex_iterator vi, vi_end, next;
	boost::tie(vi, vi_end) = boost::vertices(dual);
	for (next = vi; vi != vi_end; vi = next) {
		++next;
		if (get(name,*vi) == "SOURCE")
		{
			source = *vi;
		}
		else if (get(name,*vi) == "TARGET")
		{
			target = *vi;
		}
	}
	timeval start2, end2;
	gettimeofday(&start, 0);

	//now execute a max flow algorithm from boost lib on the dual
	//there are 3 options: edmonds-karps, push-relabel or boykov-kolmogorov
	//The latter has to be used here since it's the only one that accepts non-null capacities on an edge and at the same time on its reverse edge
	flow = boost::boykov_kolmogorov_max_flow(dual,source,target);
//	flow = boost::push_relabel_max_flow(dual, source, target);
	printf("flow is %f \n", flow);

	// now bend edges according to flow found by the algorithm
	// (flow  = capacity - residual capacity
	// arc-heights implied by flow values)

	applyFlow(&dual,&indexPrimalFaces,&indexDualVertices);

	gettimeofday(&end, 0);
	cout << start.tv_sec << ':' << start.tv_usec << endl;
	cout << end.tv_sec << ':' << end.tv_usec << endl;
	cout << "duration flow algo and applying flow: " << end.tv_sec-start.tv_sec << ":" << end.tv_usec-start.tv_usec << endl;
	store(dual,&indexPrimalFaces,&indexDualVertices,"ouput.ipe");
	printf("quit carto balance \n");

}


void MyWindow::printDualToFile(dualGraph dual)
{


	ofstream outFile("DUAL");

	property_map<dualGraph, vertex_name_t>::type name = get(vertex_name_t(), dual);
	property_map<dualGraph, edge_capacity_t>::type capacity = get(edge_capacity_t(), dual);
	property_map<dualGraph, vertex_initArea_t>::type initArea = get(vertex_initArea_t(), dual);
	property_map<dualGraph, vertex_conti_t>::type conti = get(vertex_conti_t(), dual);
	property_map<dualGraph, edge_residual_capacity_t>::type res_capacity = get(edge_residual_capacity_t(), dual);
	property_map<dualGraph, edge_reverse_t>::type reverse = get(edge_reverse_t(), dual);
	property_map<dualGraph, vertex_weight_t>::type weight = get(vertex_weight_t(), dual);

	graph_traits<dualGraph>::vertex_iterator vi, vi_end, next;
	tie(vi, vi_end) = vertices(dual);
	for (next = vi; vi != vi_end; vi = next)
	{
	    ++next;
	    //fprintf(pFile, "\n Vertex %s has weight %f \n",get(name,*vi), get(weight,*vi));

	    outFile << "\n Vertex " << get(name,*vi) << "(Continent " << get(conti,*vi) << ") has weight " << get(weight,*vi) << " and init area " << get(initArea,*vi)<< "\n";
	    graph_traits<dualGraph>::out_edge_iterator ei, ei_end, eNext;
	    tie(ei, ei_end) = out_edges(*vi, dual);
	    for (eNext = ei; ei != ei_end; ei = eNext)
		{
		//	fprintf(pFile,"edge to %s with cap %f and edge from %s with cap %f \n", get(name,target(*ei,dual)),get(capacity,*ei),get(name,target(*ei,dual)), get(capacity,get(reverse,*ei)));
			outFile << "edge to " << get(name,target(*ei,dual)) << " with cap " << get(capacity,*ei) << " and edge from " << get(name,target(*ei,dual)) << " with cap " << get(capacity,get(reverse,*ei)) << "\n";
	    	++eNext;
		}
	}

}


/* Bends the edges of the subdivision according to the flow values implied by the capacity and residual capacity map
 *
 */
void MyWindow::applyFlow(dualGraph *dual,primalMap *indPrimalFaces, dualMap *indDualVertices)
{

	primalMap indexPrimalFaces = *indPrimalFaces;

	double flow,flowSum;
	double area;
	double s;
	double h;
	int i;
	Conic_halfedge_handle hhan;

	ofstream outFile("FLOW");

	//get property maps
	property_map<dualGraph, vertex_name_t>::type name = get(vertex_name_t(), *dual);
	property_map<dualGraph, edge_capacity_t>::type capacity = get(edge_capacity_t(), *dual);
	property_map<dualGraph, edge_residual_capacity_t>::type res_capacity = get(edge_residual_capacity_t(), *dual);
	property_map<dualGraph, edge_reverse_t>::type reverse = get(edge_reverse_t(), *dual);
	property_map<dualGraph, vertex_weight_t>::type weight = get(vertex_weight_t(), *dual);
	property_map<dualGraph, vertex_initArea_t>::type initArea = get(vertex_initArea_t(), *dual);
	property_map<dualGraph, edge_primalEdges_t>::type primals = get(edge_primalEdges_t(), *dual);
	property_map<dualGraph, edge_primalEdgesCap_t>::type primalCaps = get(edge_primalEdgesCap_t(), *dual);
	property_map<dualGraph, edge_primalEdgesFlow_t>::type primalFlow = get(edge_primalEdgesFlow_t(), *dual);

	boost::graph_traits<dualGraph>::vertex_iterator u_iter, u_end;
	boost::graph_traits<dualGraph>::out_edge_iterator ei, e_end;


	 for (tie(u_iter, u_end) = vertices(*dual); u_iter != u_end; ++u_iter)
	 {
		flowSum = 0;
		for (tie(ei, e_end) = out_edges(*u_iter, *dual); ei != e_end; ++ei)
		{
				//obtain all primal edges lying on the border between the two faces
				vector<primal_edge> border = get(primals,*ei);
				vector<double> primalCapacities = get(primalCaps,*ei);


			/*	outFile << "from " <<  get(name,*u_iter) << " to " << get(name,boost::target(*ei, *dual))
						<< "\t flow: " << (get(capacity,*ei) - get(res_capacity,*ei))
						<< "\t cap: " << get(capacity,*ei)
						<< "\t edges: " << border.size() << std::endl;
*/
				flow = (get(capacity,*ei) - get(res_capacity,*ei));
				if (get(name,*u_iter)!="TARGET" && get(name,*u_iter)!="SOURCE" && get(name,boost::target(*ei, *dual))!="TARGET" && get(name,boost::target(*ei, *dual))!="SOURCE" )
				{
					flowSum += flow;

					//now increase flow over border by considering edge after edge
					for (i= 0; i < border.size();  i++)
					{

						if (flow > 0)
						{
							hhan = border[i];

							if (flow > primalCapacities[i])
							{
								area = primalCapacities[i];
							}
							else
							{
								area = flow;
							}
							get(primalFlow,*ei).push_back(area);

							h = areaToArcHeight(area,hhan->source()->point(),hhan->target()->point());
							if (debug)
							{
								std::cout << "bending edge " << hhan->source()->point() << " - " << hhan->target()->point() << " to height " << h  << std::endl;
							}
						//	std::cout << "bending edge " << hhan->source()->point() << " - " << hhan->target()->point() << " to height " << h  << std::endl;

							bendEdge(hhan->twin(),h);

							flow = flow - area;
						}
						else
						{
							get(primalFlow,*ei).push_back(0);
						}
					}
				}
		}
		outFile  	<< "Cumulated flow/area increase leaving " << get(name,*u_iter) << ": " << flowSum << std::endl
					<< "Target area increase was: " << get(weight,*u_iter)- get(initArea,*u_iter) << std::endl
					<< "Target area: " <<  get(weight,*u_iter) << std::endl
					<< "Result area: " << get(initArea,*u_iter) + flowSum << std::endl<< std::endl;
	 }


}


/* writes the arrangement (only straight edges) into an ipe file
 *
 */
void MyWindow::storeStraight(QString filename)
{
	printf("enter write straight arr to ipe \n");

	Qt_widget_demo_tab<Conic_tab_traits> *panel = static_cast<Qt_widget_demo_tab<Conic_tab_traits> *> (myBar->currentPage());

	double x1,y1,x2,y2;

	ofstream outFile(filename);

	outFile << "<?xml version=\"1.0\"?> \n";
	outFile << "<!DOCTYPE ipe SYSTEM \"ipe.dtd\"> \n";
	outFile << "<ipe version=\"70010\" creator=\"Ipe 7.0.10\"> \n";
	outFile << "<info created=\"D:20110803231932\" modified=\"D:20110803231932\"/> \n";

	//define styles (colors, pen thickness etc)
	outFile <<  "<ipestyle name=\"basic\">\n<pen name=\"light\" value=\"0.4\"/>\n<pen name=\"heavier\" value=\"0.8\"/>\n<pen name=\"fat\" value=\"1.2\"/>\n";
	outFile <<  "<pen name=\"ultrafat\" value=\"2\"/>\n<color name=\"red\" value=\"1 0 0\"/>\n<color name=\"green\" value=\"0 1 0\"/>\n";
	outFile <<  "<color name=\"blue\" value=\"0 0 1\"/>\n<color name=\"yellow\" value=\"1 1 0\"/>\n<color name=\"orange\" value=\"1 0.647 0\"/>\n";
	outFile <<  "<color name=\"gold\" value=\"1 0.843 0\"/>\n<color name=\"purple\" value=\"0.627 0.125 0.941\"/>\n<color name=\"gray\" value=\"0.745\"/>\n";
	outFile <<  "<color name=\"lightblue\" value=\"0.678 0.847 0.902\"/>\n<color name=\"lightcyan\" value=\"0.878 1 1\"/>\n";
	outFile <<  "<color name=\"lightgray\" value=\"0.827\"/>\n<color name=\"lightgreen\" value=\"0.565 0.933 0.565\"/>\n";
	outFile <<  "<color name=\"lightyellow\" value=\"1 1 0.878\"/>\n<dashstyle name=\"dashed\" value=\"[4] 0\"/>\n<dashstyle name=\"dotted\" value=\"[1 3] 0\"/>\n";
	outFile <<  "<dashstyle name=\"dash dotted\" value=\"[4 2 1 2] 0\"/>\n<dashstyle name=\"dash dot dotted\" value=\"[4 2 1 2 1 2] 0\"/>\n";
	outFile <<  "</ipestyle>\n";

	//now add paths (segments and circular arcs)
	outFile << "<page>\n<layer name=\"alpha\"/>\n<layer name=\"beta\"/>\n<view layers=\"beta alpha\" active=\"beta\"/>\n";

	Conic_halfedge_iterator hit;
	Conic_halfedge_handle hithand;
	//write all edges that border the outer face to the global vector of this class
	for (  hit = panel->m_curves_arr->halfedges_begin();  hit != panel->m_curves_arr->halfedges_end(); hit++ )
	{

		hithand = hit;
		if (hithand->face()->is_unbounded() || hithand->twin()->face()->is_unbounded())
		{
			x1 = CGAL::to_double(hithand->source()->point().x());
			y1 = CGAL::to_double(hithand->source()->point().y());
			x2 = CGAL::to_double(hithand->target()->point().x());
			y2 = CGAL::to_double(hithand->target()->point().y());
			outFile << "<path layer=\"alpha\" stroke=\"gray\" pen=\"light\">\n";
			outFile << x1 << " " << y1 << " m \n";
			outFile << x2 << " " << y2 << " l \n";
			outFile << "</path> \n";
		}
		else
		{
			x1 = CGAL::to_double(hithand->source()->point().x());
			y1 = CGAL::to_double(hithand->source()->point().y());
			x2 = CGAL::to_double(hithand->target()->point().x());
			y2 = CGAL::to_double(hithand->target()->point().y());
			outFile << "<path layer=\"alpha\" stroke=\"gray\" pen=\"light\">\n";
			outFile << x1 << " " << y1 << " m \n";
			outFile << x2 << " " << y2 << " l \n";
			outFile << "</path> \n";
		}
	}

	outFile << "</page> \n </ipe>";
	printf("quit write straight arr to svg \n");

}


/* writes the arrangement with bending into ipe file
 *
 */
void MyWindow::store(dualGraph dual,primalMap *indPrimalFaces, dualMap *indDualVertices,QString filename)
{
	Qt_widget_base_tab    *w_demo_p1 = static_cast<Qt_widget_base_tab *> (myBar->currentPage());
	primalMap indexPrimalFaces = *indPrimalFaces;

	double flow,rflow;
	double h;
	double x1,y1,x2,y2,x3,y3,xc,yc,s;
	double diffX,diffY,diffXc,diffYc,alpha,radius,bigA,bigB,bigC;
	int i;

	Conic_halfedge_handle hhan;
	Arr_conic_point_2 source,target;

	ofstream outFile(filename);

/*
	outFile << "<?xml version=\"1.0\" standalone=\"no\"?> \n";
	outFile << "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\" \n \"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\"> \n";
	outFile << "<svg width=\"cm\" height=\"4cm\" viewBox=\" " <<  w_demo_p1->bbox.xmin() << " " <<  w_demo_p1->bbox.ymin()
			<< " " << w_demo_p1->bbox.xmax()-w_demo_p1->bbox.xmin() << " " << w_demo_p1->bbox.ymax()-w_demo_p1->bbox.ymin() << "\""
			<< " \n xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\"> \n ";
	outFile << "<title>Transformed Map</title> \n";
	outFile << "<desc>The resulting Circular-arc Cartogram </desc> \n";
*/
	outFile << "<?xml version=\"1.0\"?> \n";
	outFile << "<!DOCTYPE ipe SYSTEM \"ipe.dtd\"> \n";
	outFile << "<ipe version=\"70010\" creator=\"Ipe 7.0.10\"> \n";
	outFile << "<info created=\"D:20110803231932\" modified=\"D:20110803231932\"/> \n";

	//define styles (colors, pen thickness etc)
	outFile <<  "<ipestyle name=\"basic\">\n<pen name=\"heavier\" value=\"0.8\"/>\n<pen name=\"fat\" value=\"1.2\"/>\n";
	outFile <<  "<pen name=\"ultrafat\" value=\"2\"/>\n<color name=\"red\" value=\"1 0 0\"/>\n<color name=\"green\" value=\"0 1 0\"/>\n";
	outFile <<  "<color name=\"blue\" value=\"0 0 1\"/>\n<color name=\"yellow\" value=\"1 1 0\"/>\n<color name=\"orange\" value=\"1 0.647 0\"/>\n";
	outFile <<  "<color name=\"gold\" value=\"1 0.843 0\"/>\n<color name=\"purple\" value=\"0.627 0.125 0.941\"/>\n<color name=\"gray\" value=\"0.745\"/>\n";
	outFile <<  "<color name=\"lightblue\" value=\"0.678 0.847 0.902\"/>\n<color name=\"lightcyan\" value=\"0.878 1 1\"/>\n";
	outFile <<  "<color name=\"lightgray\" value=\"0.827\"/>\n<color name=\"lightgreen\" value=\"0.565 0.933 0.565\"/>\n";
	outFile <<  "<color name=\"lightyellow\" value=\"1 1 0.878\"/>\n<dashstyle name=\"dashed\" value=\"[4] 0\"/>\n<dashstyle name=\"dotted\" value=\"[1 3] 0\"/>\n";
	outFile <<  "<dashstyle name=\"dash dotted\" value=\"[4 2 1 2] 0\"/>\n<dashstyle name=\"dash dot dotted\" value=\"[4 2 1 2 1 2] 0\"/>\n";
	outFile <<  "</ipestyle>\n";

	//now add paths (segments and circular arcs)
	outFile << "<page>\n<layer name=\"alpha\"/>\n<layer name=\"beta\"/>\n<view layers=\"beta alpha\" active=\"beta\"/>\n";

	//get property maps
	property_map<dualGraph, vertex_name_t>::type name = get(vertex_name_t(), dual);
	property_map<dualGraph, edge_capacity_t>::type capacity = get(edge_capacity_t(), dual);
	property_map<dualGraph, edge_residual_capacity_t>::type res_capacity = get(edge_residual_capacity_t(), dual);
	property_map<dualGraph, edge_reverse_t>::type reverse = get(edge_reverse_t(), dual);
	property_map<dualGraph, vertex_weight_t>::type weight = get(vertex_weight_t(), dual);
	property_map<dualGraph, vertex_initArea_t>::type initArea = get(vertex_initArea_t(), dual);
	property_map<dualGraph, edge_primalEdges_t>::type primals = get(edge_primalEdges_t(), dual);
	property_map<dualGraph, edge_primalEdgesCap_t>::type primalCaps = get(edge_primalEdgesCap_t(), dual);
	property_map<dualGraph, edge_primalEdgesFlow_t>::type primalFlow = get(edge_primalEdgesFlow_t(), dual);
	property_map<dualGraph, edge_primalEdgesSources_t>::type primalSources = get(edge_primalEdgesSources_t(), dual);
	property_map<dualGraph, edge_primalEdgesTargets_t>::type primalTargets = get(edge_primalEdgesTargets_t(), dual);

	 boost::graph_traits<dualGraph>::vertex_iterator u_iter, u_end;
	 boost::graph_traits<dualGraph>::out_edge_iterator ei, e_end;

	 for (tie(u_iter, u_end) = vertices(dual); u_iter != u_end; ++u_iter)
	 {
		//std::cout << "start vertex " << get(name,*u_iter) << std::endl;
		for (tie(ei, e_end) = out_edges(*u_iter, dual); ei != e_end; ++ei)
		{
		//	 std::cout << "goal vertex " << get(name,boost::target(*ei, dual)) << std::endl;

			//obtain all primal edges lying on the border between the two faces
			vector<primal_edge> border = get(primals,*ei);
			vector<Arr_conic_point_2> sources = get(primalSources,*ei);
			vector<Arr_conic_point_2> targets = get(primalTargets,*ei);
			vector<double> primalCapacities = get(primalCaps,*ei);
			vector<double> flows = get(primalFlow,*ei);
			vector<double> reverseFlows = get(primalFlow,get(reverse,*ei));

			flow = get(capacity,*ei) - get(res_capacity,*ei);
	//		std::cout<< "flow here is " << flow << " and there are " << border.size() << " edges on this border " << std::endl;

			//now increase flow over border by considering edge after edge
			if  (get(name,*u_iter)!="TARGET" && get(name,*u_iter)!="SOURCE" && get(name,boost::target(*ei, dual))!="TARGET" && get(name,boost::target(*ei, dual))!="SOURCE")
			{
				for (i=0; i < border.size();i++)
				{
					hhan = border[i];
					//source and target need to be switched since in function apply flow, the flow is always sent over twin edge
					source = targets[i];
					target = sources[i];
					flow = flows[i];
					rflow = reverseFlows[i];
				//	std::cout<< "primal edge " << source << " to " << target << " with flow of " << flow << std::endl;

					if (flow > 0)
					{
						h = areaToArcHeight(flow,source,target);
						x1 = CGAL::to_double(source.x());
						y1 = CGAL::to_double(source.y());
						x2 = CGAL::to_double(target.x());
						y2 = CGAL::to_double(target.y());
					 	s = sqrt(pow(x1-x2,2) + pow(y1-y2,2));

						K_point sourceK (x1,y1);
						K_point targetK (x2,y2);
						K_point edgeCenter = sourceK + (targetK - sourceK) / 2.0;

						//now do perpendicular projection from edge center in the right direction (into the polygon)
						//the inside of the polygon lies on the right side of cc
						//alpha is the angle of the contour edge with the global x-axis
						alpha = atan(fabs(y2-y1)/fabs(x2-x1));
						diffX = fabs(cos(PI/2.0-alpha)*h);
						diffY = fabs(sin(PI/2.0-alpha)*h);

						//case distinction depends on in which quadrant the halfedge vector
						//lies when its source is the source of the coordinate-system
						K_point middle;

						if (x2 <= x1 && y2 <= y1)
						{
							middle = K_point(edgeCenter.x()+diffX,edgeCenter.y()-diffY);
						}
						else if (x2 >= x1 && y2 >=y1)
						{
							middle = K_point(edgeCenter.x()-diffX,edgeCenter.y()+diffY);
						}
						else if (x2 >= x1 && y2 <= y1)
						{
							middle = K_point(edgeCenter.x()+diffX,edgeCenter.y()+diffY);
						}
						else
						{
							middle = K_point(edgeCenter.x()-diffX,edgeCenter.y()-diffY);
						}
				//		std::cout << "storing the arc through 3 points " << source << " = " << sourceK << "   /   " << middle << "   /   " << target << " = " << targetK << std::endl;

						x3 = middle.x();
						y3 = middle.y();


						bigA = sqrt(pow(x1-x2,2)+pow(y1-y2,2));
						bigB = sqrt(pow(x1-x3,2)+pow(y1-y3,2));
						bigC = sqrt(pow(x2-x3,2)+pow(y2-y3,2));

						radius = bigA*bigB*bigC / sqrt((bigA+bigB+bigC)*(-bigA+bigB+bigC)*(bigA-bigB+bigC)*(bigA+bigB-bigC)) ;



						xc = ((pow(y1,2)+pow(x1,2))*(y2-y3)+(pow(y2,2)+pow(x2,2))*(y3-y1)+(pow(y3,2)+pow(x3,2))*(y1-y2)) / (2*(x1*(y2-y3)+x2*(y3-y1)+x3*(y1-y2)));
						yc = ((pow(y1,2)+pow(x1,2))*(x3-x2)+(pow(y2,2)+pow(x2,2))*(x1-x3)+(pow(y3,2)+pow(x3,2))*(x2-x1)) / (2*(x1*(y2-y3)+x2*(y3-y1)+x3*(y1-y2)));
						K_point circleCenter(xc,yc);

			//			std::cout << "gives radius " << radius << " and center " << circleCenter << std::endl;

						outFile << "<path layer=\"beta\" stroke=\"black\" pen=\"fat\">\n";
						outFile << x1 << " " << y1 << " m \n";
						if (h >= s/2)
						{
							outFile << -radius << " 0 0 " << radius << " " << xc << " " << yc << " " << x2 << " " << y2 << " a \n";
						}
						else
						{
							outFile << radius << " 0 0 " <<  -radius << " " << xc << " " << yc << " " << x2 << " " << y2 << " a \n";

						}
						outFile << "</path> \n";
						outFile << "<path layer=\"alpha\" stroke=\"gray\" pen=\"fat\">\n";
						outFile << x1 << " " << y1 << " m \n";
						outFile << x2 << " " << y2 << " l \n";
						outFile << "</path> \n";

					}
					else if (flow == 0 and rflow == 0)
					{
						x1 = CGAL::to_double(source.x());
						y1 = CGAL::to_double(source.y());
						x2 = CGAL::to_double(target.x());
						y2 = CGAL::to_double(target.y());

						//outFile << "<path d=\"M" << x1 << " " << y1 << " L " << x2 << " " <<y2 << "\" stroke=\"blue\" stroke-width=\"1\" /> \n ";
						outFile << "<path layer=\"beta\" stroke=\"black\" pen=\"fat\">\n";
						outFile << x1 << " " << y1 << " m \n";
						outFile << x2 << " " << y2 << " l \n";
						outFile << "</path> \n";
					}


				}
			}
		}
	}
	outFile << "</page> \n </ipe>";

}


//reads weight data from a file
void MyWindow::cartogram_weights()
{

	  Qt_widget_demo_tab<Conic_tab_traits> *panel =
	  		static_cast<Qt_widget_demo_tab<Conic_tab_traits> *> (myBar->currentPage());
	  Conic_face_iterator fit;
	  fit = panel->m_curves_arr->faces_begin();
	  QString filename = QFileDialog::getOpenFileName(QString::null, 0, this,"weight file open", "open Weight File" );
	  if ( !filename.isEmpty() )
	  {
		updateMode( dragMode );
	  }
	  else
	    statusBar()->message( "File Open abandoned", 2000 );
	  printf("Open Weight File \n");
	  std::ifstream inputFile(filename.ascii());
	  if (! inputFile.is_open())
	  {
	    std::cout << "Error opening weight file" << std::endl;
	    return;
	  }

	  char* value;
	  char* country;
	  string line;
	  double area;


	  while (inputFile.good() && (fit != panel->m_curves_arr->faces_end()))
	  {
		  if (!fit->is_unbounded())
		  {
			  getline(inputFile,line);
			  int firstSpace = line.find(" ");
			  value = new char [firstSpace];
			  strcpy(value,line.substr(0,firstSpace).c_str());
			  country = new char[line.length()-firstSpace];
			  strcpy(country,line.substr(firstSpace+1,line.length()).c_str());
			  area =  FaceArea(fit);
			  fit->setWeight(atof(value));
			  fit->setInitialArea(area);
			  fit->setCountry(country);
		  }
		  else
		  {
			  fit->setWeight(0);
			  fit->setInitialArea(0);
			  fit->setCountry((char *)"OUTERFACE");
		  }
		  fit++;
	  }
}



//sets up the weights for all regions such that the total
//area remains the same but each region is weighted with the given data tuple
void MyWindow::setUpWeights()
{
	Qt_widget_demo_tab<Conic_tab_traits> *panel =
				static_cast<Qt_widget_demo_tab<Conic_tab_traits> *> (myBar->currentPage());
	//dont do it unless we have more than one bounded face (otherwise we dont need to scale the weights)
	if (panel->m_curves_arr->number_of_faces() > 2)
	{
		Conic_face_iterator fit;
		Conic_face_handle fhan;
		float totalWeight = 0.0;
		float totalOriginalArea = 0.0;
		// Calculate the total original area and total target area of the subdivision
		for (fit = panel->m_curves_arr->faces_begin(); fit != panel->m_curves_arr->faces_end(); ++fit)
		{
		   fhan = fit;
		   if (!fhan->is_unbounded())
		   {
			   totalWeight += fhan->getWeight();
			   totalOriginalArea += FaceArea(fhan);
		   }
		}
		// set new target areas weighted by the quotient of total original area and total target area (we preserve the original area)
		for (fit = panel->m_curves_arr->faces_begin(); fit != panel->m_curves_arr->faces_end(); ++fit)
		{
		   fhan = fit;
		   if (!fhan->is_unbounded())
			   fhan->setWeight(totalOriginalArea*((fhan->getWeight())/totalWeight));
			if (debug)
			{
				std::cout << fhan->getCountry() << " has weight " << fhan->getWeight() << " and init area " << fhan->getInitialArea() << std::endl;
			}

		}
	}

}

//calculates the area of the face
//if the face has non-colinear edges (e.g. conics) then an approximation with precision given as parameter will
//be  done
float MyWindow::FaceArea(Conic_face_handle f)
{
	Qt_widget_demo_tab<Conic_tab_traits> *panel =
			static_cast<Qt_widget_demo_tab<Conic_tab_traits> *> (myBar->currentPage());
    Conic_ccb_halfedge_circulator cc=f->outer_ccb();
	bool colinearPolygon = true;
	do	{
	   colinearPolygon = colinearPolygon && (cc->curve().orientation()==0);
	} while (++cc != f->outer_ccb());

	if (colinearPolygon)
	{
		return extractRoughPolygonFromFace(f).area();
	}
	else
	{
		return extractFinePolygonFromFace(f,1).area();
	}

}



//returns the polygon representation from a given face handle
//does simply take endpoints no matter what type the segment has (conic, polyline..)
My_polygon MyWindow::extractRoughPolygonFromFace(Conic_face_handle f)
{
	Qt_widget_demo_tab<Conic_tab_traits> *panel =
		static_cast<Qt_widget_demo_tab<Conic_tab_traits> *> (myBar->currentPage());
    std::list< Coord_point > pts; // Polygon points
    Conic_ccb_halfedge_circulator cc=f->outer_ccb();
    do {
      if (panel->antenna(cc))
        continue;

      Conic_halfedge_handle hhan = cc;
      // Get the coordinates of the curve's source and target
      double sx = CGAL::to_double(hhan->source()->point().x()),
             sy = CGAL::to_double(hhan->source()->point().y()),
             tx = CGAL::to_double(hhan->target()->point().x()),
             ty = CGAL::to_double(hhan->target()->point().y());

      Coord_point coord_source(sx,sy);
      Coord_point coord_target(tx,ty);

      pts.push_back(coord_source);

    } while (++cc != f->outer_ccb());

    // make polygon from the outer ccb of the face f
    My_polygon pgn (pts.begin() , pts.end());

    return pgn;
}



//returns the polygon representation from a given face handle
//does approximate non colinear curves with a given precision
My_polygon MyWindow::extractFinePolygonFromFace(Conic_face_handle f, float precision)
{
	Qt_widget_demo_tab<Conic_tab_traits> *panel =
		static_cast<Qt_widget_demo_tab<Conic_tab_traits> *> (myBar->currentPage());
    std::list< Coord_point > pts; // Polygon points
    Conic_ccb_halfedge_circulator cc=f->outer_ccb();
    do {
      if (panel->antenna(cc))
        continue;

      Conic_halfedge_handle hhan = cc;
      Arr_xconic_2 c = hhan->curve();
      // Get the coordinates of the curve's source and target
      double sx = CGAL::to_double(hhan->source()->point().x()),
             sy = CGAL::to_double(hhan->source()->point().y()),
             tx = CGAL::to_double(hhan->target()->point().x()),
             ty = CGAL::to_double(hhan->target()->point().y());

      Coord_point coord_source(sx,sy);
      Coord_point coord_target(tx,ty);

      if (c.orientation() == CGAL::COLLINEAR)
		 pts.push_back(coord_source );
      else
      {
		   // If the curve is monotone, than its source and its target have the
		   // extreme x coordinates
		   bool is_source_left = (sx < tx);
		   int  x_min = is_source_left ? panel->x_pixel(sx) : panel->x_pixel(tx);
		   int  x_max = is_source_left ? panel->x_pixel(tx) : panel->x_pixel(sx);
		   double curr_x, curr_y;
		   int  x;
		   Arr_conic_point_2 px;

		   pts.push_back(coord_source );

		   if (is_source_left) {
			 for (x = x_min + precision; x < x_max; x+=precision) {
			   curr_x = panel->x_real(x);
			   Alg_kernel   ker;
			   Arr_conic_point_2 curr_p(curr_x, 0);
			   if (!(ker.compare_x_2_object()(curr_p, c.left()) != CGAL::SMALLER &&
					   ker.compare_x_2_object()(curr_p, c.right()) != CGAL::LARGER))
				   continue;
			   px = c.point_at_x (curr_p);
			   curr_y = CGAL::to_double(px.y());
			   pts.push_back(Coord_point(curr_x ,curr_y));
			 }
		   }
		   else {
			 for (x = x_max; x > x_min; x-=precision) {
			   curr_x = panel->x_real(x);
			   Alg_kernel   ker;
			   Arr_conic_point_2 curr_p(curr_x, 0);
			   if (!(ker.compare_x_2_object() (curr_p, c.left()) != CGAL::SMALLER &&
					 ker.compare_x_2_object() (curr_p, c.right()) != CGAL::LARGER))
				   continue;
			   px = c.point_at_x (Arr_conic_point_2(curr_x, 0));
			   curr_y = CGAL::to_double(px.y());
			   pts.push_back(Coord_point(curr_x,curr_y));
			 }
		   }
		   pts.push_back(coord_target );
		}
	   } while (++cc != f->outer_ccb());

    // make polygon from the outer ccb of the face f
    My_polygon pgn (pts.begin() , pts.end());

    return pgn;
}




void MyWindow::cartogram_it()
{
}


// traverses all bounded faces and pops up their weights in an extra window
// erase later, won't be needed
void MyWindow::print_all_weights()
{
	 QColor cRed = Qt::red;
	QColor cBlack = Qt::black;
	Qt_widget_demo_tab<Conic_tab_traits> *panel =static_cast<Qt_widget_demo_tab<Conic_tab_traits> *> (myBar->currentPage());
	Conic_face_handle ubf = panel->m_curves_arr->unbounded_face();
	std::cout << panel->m_curves_arr->number_of_faces() << " faces:" << std::endl;
	Conic_face_iterator fit;
	Conic_face_handle fhan;
	for (fit = panel->m_curves_arr->faces_begin(); fit != panel->m_curves_arr->faces_end(); ++fit)
	{
			fhan = fit;  //*(dynamic_cast<Conic_face_handle*>(&*fit));
			if (!fhan->is_unbounded())
			{
					panel->set_face_color(fhan, cRed);
					something_changed();
					Notifier *displayWeight = new Notifier(fhan->getWeight(), myBar , this ,number_of_tabs ,panel , m_scailing_factor , colors_flag);
					displayWeight->move (0, 0);
					if ( displayWeight->exec() )
					{
							int i=0;
					}
					delete displayWeight;
					float area =  fhan->getWeight();
					printf("area: ");
					printf("%f \n",area);
					printf("\n");
					panel->set_face_color(fhan, cBlack);
					something_changed();
			}
	}


}

//same as 'extractRoughPolygonFromFace' but different return type (other geometry kernel for polygon)
//look for a maybe more elegant way of doing this (handling template in return type instead of 2 equivalent functions)
My_polygon_K MyWindow::extractPolygonKFromFace(Conic_face_handle f)
{
	Qt_widget_demo_tab<Conic_tab_traits> *panel =
		static_cast<Qt_widget_demo_tab<Conic_tab_traits> *> (myBar->currentPage());
    std::list< K_point > pts; // Polygon points
    Conic_ccb_halfedge_circulator cc=f->outer_ccb();
    do {
      if (panel->antenna(cc))
        continue;

      Conic_halfedge_handle hhan = cc;
      // Get the coordinates of the curve's source and target
      double sx = CGAL::to_double(hhan->source()->point().x()),
             sy = CGAL::to_double(hhan->source()->point().y());

      K_point coord_source(sx,sy);

      pts.push_back(coord_source);

    } while (++cc != f->outer_ccb());

    // make polygon from the outer ccb of the face f
	My_polygon_K pgn (pts.begin() , pts.end());
	return pgn;
}




//Flips the straight edges corresponding to the halfedge hhan so that the circular arc has height h
//the edge will be bent to the right side of hhan (if hhan is a halfedge on the outer boundary of a polygon
//then the edge will be flipped inwards w.r.t this polygon)
void MyWindow::bendEdge(Conic_halfedge_handle hhan, double h)
{

	if (h > 0 )
	{
		Qt_widget_demo_tab<Conic_tab_traits> *panel =
					static_cast<Qt_widget_demo_tab<Conic_tab_traits> *> (myBar->currentPage());


		// if using modify_edge then the arcs has to be x-monotone (first line), otherwise we need the second line
		//Arr_xconic_2 newCurve  (static_cast<Arr_conic_2>(obtainCircularArc(hhan,h)));
		Arr_conic_2 newCurve(obtainCircularArc(hhan,h));


		// if no non x-monotone curves occur, we can use the set_curve function (one of the two lines below)
		//hhan->twin()->set_curve(newCurve);
		//panel->m_curves_arr->modify_edge(hhan,newCurve);

		// if non x-monotone curves occur, we cannot use the set_curve function -> we need to remove and reinsert the halfedges with the new curve
		panel->m_curves_arr->remove_edge(hhan);
		insert(*(panel->m_curves_arr), newCurve);

		something_changed();
	}

}


/* creates circular arc under straight edge hhan with height h
 *
 */
Arr_conic_2 MyWindow::obtainCircularArc(Conic_halfedge_handle hhan, double h)
{
	//first store source and target point of edge (need rationals for arc construction)
	//Strange error when creating rationals: always 1/10000 less then float given (occurs when casting statically float to integer)
	//provisional fix: just add 0.5 to float value in order to avoid incorrect rounding (probably the source of imprecision when casting)
	Rational x1 (static_cast<int> (10000*CGAL::to_double(hhan->source()->point().x())+0.5), 10000);
	Rational x2 (static_cast<int> (10000*CGAL::to_double(hhan->target()->point().x())+0.5), 10000);
	Rational y1 (static_cast<int> (10000*CGAL::to_double(hhan->source()->point().y())+0.5), 10000);
	Rational y2 (static_cast<int> (10000*CGAL::to_double(hhan->target()->point().y())+0.5), 10000);
	Rat_point_2 source (x1,y1);
	Rat_point_2 target (x2,y2);


 /* debugging:
  * printf("old point p1 was %f / %f \n", CGAL::to_double(hhan->source()->point().x()), CGAL::to_double(hhan->source()->point().y()));
	printf("new point p1  is %f / %f \n",  CGAL::to_double(x1),CGAL::to_double(y1));
	printf("old point p2 was %f / %f \n", CGAL::to_double(hhan->target()->point().x()), CGAL::to_double(hhan->target()->point().y()));
	printf("new point p2  is %f / %f \n",  CGAL::to_double(x2),CGAL::to_double(y2));
  */
	//derive 3rd point of new circular arc (perpendicular projection of middle point of edge by h into the polygon)
	//first get the middle point of the edge
	Rat_point_2 edgeCenter = source + (target - source) / 2.0;

	//now do perpendicular projection from edge center in the right direction (into the polygon)
	//the inside of the polygon lies on the right side of cc
	//alpha is the angle of the contour edge with the global x-axis
	double alpha = atan(fabs(CGAL::to_double(y2)-CGAL::to_double(y1))/fabs(CGAL::to_double(x2)-CGAL::to_double(x1)));
	Rational diffX  (static_cast<int> (fabs(10000*cos(PI/2.0-alpha)*h)),10000);
	Rational diffY  (static_cast<int> (fabs(10000*sin(PI/2.0-alpha)*h)),10000);

	//case distinction depends on in which quadrant the halfedge vector
	//lies when its source is the source of the coordinate-system
	Rat_point_2 middle;
	if (hhan->twin()->source()->point().x() <= hhan->twin()->target()->point().x() && hhan->twin()->source()->point().y() <= hhan->twin()->target()->point().y())
	{
		middle = Rat_point_2(edgeCenter.x()+diffX,edgeCenter.y()-diffY);
	}
	else if (hhan->twin()->source()->point().x() >= hhan->twin()->target()->point().x() && hhan->twin()->source()->point().y() >= hhan->twin()->target()->point().y())
	{
		middle = Rat_point_2(edgeCenter.x()-diffX,edgeCenter.y()+diffY);
	}
	else if (hhan->twin()->source()->point().x() >= hhan->twin()->target()->point().x() && hhan->twin()->source()->point().y() <= hhan->twin()->target()->point().y())
	{
		middle = Rat_point_2(edgeCenter.x()+diffX,edgeCenter.y()+diffY);
	}
	else
	{
		middle = Rat_point_2(edgeCenter.x()-diffX,edgeCenter.y()-diffY);
	}

/*  // debugging: show the vertices obtained from above calculation
	Arr_conic_point_2 cen (CGAL::to_double(edgeCenter.x()),CGAL::to_double(edgeCenter.y()));
	panel->m_curves_arr->insert_in_face_interior(cen,fit);
	Arr_conic_point_2 mid (CGAL::to_double(middle.x()),CGAL::to_double(middle.y()));
	panel->m_curves_arr->insert_in_face_interior(mid,fit);
*/

	return Arr_conic_2 (source,middle,target);
}


/*
 * displays all edges of the straight skeleton for the face which is passed to this function
 */
void MyWindow::showSkeleton(Conic_face_handle fhan)
{
	printf("enter show skeleton \n");
	Qt_widget_demo_tab<Conic_tab_traits> *panel = static_cast<Qt_widget_demo_tab<Conic_tab_traits> *> (myBar->currentPage());
	Rat_segment_2 segment;
	Halfedge_iterator i;
	SsPtr skeleton;

	if (!fhan->is_unbounded())
	{
		My_polygon_K pgn = extractPolygonKFromFace(fhan);
		skeleton = CGAL::create_interior_straight_skeleton_2(pgn.vertices_begin(), pgn.vertices_end());

	}
	else
	{
		printf("unbounded \n");

		My_polygon_K outer = getBoundingPolygon();
		My_polygon_with_holes_K poly(outer);
		for (int i = 0; i < numberOfContinents; i++)
		{
			printf("iter \n");
			poly.add_hole(getContinentPolygon(i));

		}
		skeleton = CGAL::create_interior_straight_skeleton_2(poly);
	}

	printf("adding edges \n");

	for (  i = (*skeleton).halfedges_begin(); i != (*skeleton).halfedges_end(); ++i )
	{
		//printf("Adding straight segment %f / %f  to %f / %f \n",CGAL::to_double(i->vertex()->point().x()),CGAL::to_double(i->vertex()->point().y()),CGAL::to_double(i->opposite()->vertex()->point().x()),CGAL::to_double(i->opposite()->vertex()->point().y()));
		segment = Rat_segment_2(Rat_point_2(i->vertex()->point().x(),i->vertex()->point().y()),Rat_point_2(i->opposite()->vertex()->point().x(),i->opposite()->vertex()->point().y()));
		insert(*(panel->m_curves_arr),Arr_conic_2(segment));
	}
	something_changed();
}

void MyWindow::printSkeleton(SsPtr skeleton)
{
	printf("enter print skeleton \n");
	Qt_widget_demo_tab<Conic_tab_traits> *panel = static_cast<Qt_widget_demo_tab<Conic_tab_traits> *> (myBar->currentPage());
	Rat_segment_2 segment;
	Halfedge_iterator i;

	for (  i = (*skeleton).halfedges_begin(); i != (*skeleton).halfedges_end(); ++i )
	{
		segment = Rat_segment_2(Rat_point_2(i->vertex()->point().x(),i->vertex()->point().y()),Rat_point_2(i->opposite()->vertex()->point().x(),i->opposite()->vertex()->point().y()));
		insert(*(panel->m_curves_arr),Arr_conic_2(segment));
	}
	something_changed();
	printf("quit print skeleton \n");

}


/*
 * displays edges of straight skeleton face that is associated with contour edge passed as parameter to this function
 */
void MyWindow::showSkeletonFace(Conic_face_handle fhan,Conic_halfedge_handle hhan)
{
	Qt_widget_demo_tab<Conic_tab_traits> *panel =
					static_cast<Qt_widget_demo_tab<Conic_tab_traits> *> (myBar->currentPage());
	Rat_segment_2 segment;
	Halfedge_iterator i;
	bool found = false;
	My_polygon_K pgn = extractPolygonKFromFace(fhan);
	SsPtr skeleton = CGAL::create_interior_straight_skeleton_2(pgn.vertices_begin(), pgn.vertices_end());
	// retrieve the skeleton edge which corresponds to the original edge we are looking for
	for (  i = (*skeleton).halfedges_begin(); (found == false) && (i != (*skeleton).halfedges_end()); ++i )
	{
		if (!(i->is_bisector()) && compareEdges(hhan,i))
		{	for (  i = (*skeleton).halfedges_begin(); i != (*skeleton).halfedges_end(); ++i )

		   found = true;
		}
	}
	//now traverse all skeleton edges that border the face belonging to the skeleton contour edge hhan (other than actual border contour edge)
	Halfedge_handle skeletonFaceIterator;
	for ( skeletonFaceIterator = i->opposite()->face()->halfedge()->next(); skeletonFaceIterator != i->opposite()->face()->halfedge(); skeletonFaceIterator = skeletonFaceIterator->next())
	{
		segment = Rat_segment_2(Rat_point_2(skeletonFaceIterator->vertex()->point().x(),skeletonFaceIterator->vertex()->point().y()),Rat_point_2(skeletonFaceIterator->opposite()->vertex()->point().x(),skeletonFaceIterator->opposite()->vertex()->point().y()));
		insert(*(panel->m_curves_arr),Arr_conic_2(segment));
	}
	something_changed();
}



//returns the estimated maximum area change that can be done by bending the edge corresponding to hhan inside the corresponding skeleton region
//which borders the face fhan inwards (circular arc maximally bent in straight skeleton region)
//precondition: hhan must be part of fhan and hhan is a straight line segment (coolinear is true)
//this is only an APPROXIMATIOn (smaller angle of thwo neighboring skeleton face edges)
double MyWindow::calcBendingCapacityEstimate(Conic_face_handle fhan,Conic_halfedge_handle hhan )
{
	Qt_widget_demo_tab<Conic_tab_traits> *panel =
				static_cast<Qt_widget_demo_tab<Conic_tab_traits> *> (myBar->currentPage());

	double h;
	double A;

	if (!fhan->is_unbounded())
	{
		My_polygon_K pgn = extractPolygonKFromFace(fhan);
		SsPtr straightSkeleton;
		bool found=false;
		double h;
		// calculate the straight skeleton of the face if this has not been done previously
		straightSkeleton = CGAL::create_interior_straight_skeleton_2(pgn.vertices_begin(), pgn.vertices_end());
		Halfedge_iterator i;
		for (i = (*straightSkeleton).halfedges_begin(); (found == false) && (i != (*straightSkeleton).halfedges_end()); ++i)
		{
			if (!(i->is_bisector()) && compareEdges(hhan,i))
			{
			   //get the two relevant bisectors of skeleton face
			   Halfedge_handle bisector1 = i->vertex()->primary_bisector();
			   Halfedge_handle bisector2 = i->opposite()->vertex()->primary_bisector();
			   //calculate the two angles between contour edge i and the bisectors
			   double alpha = calculateAngle(bisector1, i, bisector1->opposite()->vertex()->point());
			   double beta = calculateAngle(bisector2, i,bisector2->opposite()->vertex()->point());
			   //get smaller angle (the one which limits bending)
			   if (alpha > beta)
			   {
				   alpha = beta;
			   }
			   //calculate length of edge hhan (is a straight line)
			   double s = sqrt(fabs(CGAL::to_double(hhan->source()->point().x()) - CGAL::to_double(hhan->target()->point().x()))
							   *fabs(CGAL::to_double(hhan->source()->point().x()) - CGAL::to_double(hhan->target()->point().x()))
							   +fabs(CGAL::to_double(hhan->source()->point().y()) - CGAL::to_double(hhan->target()->point().y()))
							   *fabs(CGAL::to_double(hhan->source()->point().y()) - CGAL::to_double(hhan->target()->point().y())));
			   //apply formula for circular segment area (base angle theta = 2*alpha)
			   h = s/2*tan(alpha/2); //height
			   double r = (4*h*h+s*s)/(8*h); //radius
			   A = r*r/2*(2*alpha-sin(2*alpha)); //area
			   //debugging
			   //printf("alpha: %f \n s: %f \n h: %f \n r: %f \n A: %f \n",alpha,s,h,r,A);
			   //end this loop with this iteration (because correct contour edge has been found)
			   found = true;
			}
		}
	}
	else
	{

	   //calculate the two angles between contour edge i and the two neighbors on outer face
	   double alpha = PI/4;
	   //calculate length of edge hhan (is a straight line)
	   double s = sqrt(fabs(CGAL::to_double(hhan->source()->point().x()) - CGAL::to_double(hhan->target()->point().x()))
					   *fabs(CGAL::to_double(hhan->source()->point().x()) - CGAL::to_double(hhan->target()->point().x()))
					   +fabs(CGAL::to_double(hhan->source()->point().y()) - CGAL::to_double(hhan->target()->point().y()))
					   *fabs(CGAL::to_double(hhan->source()->point().y()) - CGAL::to_double(hhan->target()->point().y())));
	   //apply formula for circular segment area (base angle theta = 2*alpha)
	   h = s/2*tan(alpha/2); //height
	   double r = (4*h*h+s*s)/(8*h); //radius
	   A = r*r/2*(2*alpha-sin(2*alpha)); //area
	}

	return A;
}



//returns the exact maximum area change that can be done by bending the edge corresponding inside the corresponding skeleton region to hhan
//which borders the face fhan inwards (circular arc maximally bent in straight skeleton region)
//precondition: hhan must be part of fhan and hhan is a straight line segment (coolinear is true)
//this is the EXACT computation of the bending height
double MyWindow::calcBendingCapacityExact(Conic_face_handle fhan,Conic_halfedge_handle hhan,primalMap *indPrimalFaces )
{
		Qt_widget_demo_tab<Conic_tab_traits> *panel = static_cast<Qt_widget_demo_tab<Conic_tab_traits> *> (myBar->currentPage());

		//calculate the straight skeleton of the face if this has not been done previously
		My_polygon_K pgn;
		My_polygon_with_holes_K pgnHoles;
		SsPtr straightSkeleton;

		Skeleton_face skeletonFace;
		Halfedge_handle skeletonFaceIterator;
		bool found = false;
		double h,hprov,A,s = 0;
		double alpha,beta,m1,m2,m3,m4,x1,x2,x3,x4,y1,y2,y3,y4,xc,yc;
		double xx1,xx2,xx3,xx4,yy1,yy2,yy3,yy4;
		double xxx1,xxx2,xxx3,xxx4,yyy1,yyy2,yyy3,yyy4;
		double a,b,b1,b2,n,m = 0;
		double aux1, aux2;
		double yres1,yres2,xres1,xres2;
		double la,lb,lc;
		int i;
		Arr_conic_point_2 p1,p2,p3,p4,p5,p6;

		// p1 is source point and p2 is target ppoint of first segement
		x1 = CGAL::to_double(hhan->source()->point().x());
		y1 = CGAL::to_double(hhan->source()->point().y());
		x2 = CGAL::to_double(hhan->target()->point().x());
		y2 = CGAL::to_double(hhan->target()->point().y());
		// s is its length
		s =  sqrt(pow(x1 - x2,2)  + pow(y1-y2,2));

		// get rotation matrix such that after its application the edge is horizontal, points from left to right and the associated face lies above edge (to its left)
		if (x1 == x2)
		{
			//vertical case
			if (y1 > y2)
			{
				//case where we have to rotate clockwise
				alpha = PI/2;;
	//				printf("Rotate %f degrees counterclockwise \n",alpha/(2*PI)*360);
				m1 = cos(alpha);
				m2 = -sin(alpha);
				m3 = sin (alpha);
				m4 = cos(alpha);

			}
			else
			{
				//case where we have to rotate counterclockwise
				alpha = PI/2;
			//	printf("Rotate %f degrees clockwise \n",alpha/(2*PI)*360);
				m1 = cos(alpha);
				m2 = sin(alpha);
				m3 = -sin (alpha);
				m4 = cos(alpha);
			}
		}
		else if (x1 < x2)
		{
			//case where source.x is smaller then target.x, means that associated face lies above edge initially
			if (y1 <= y2)
			{
				//case where we have to rotate clockwise
				alpha = atan((y2-y1)/(x2-x1));
		//		printf("Rotate %f degrees clockwise \n",alpha/(2*PI)*360);
				m1 = cos(alpha);
				m2 = sin(alpha);
				m3 = -sin (alpha);
				m4 = cos(alpha);

			}
			else
			{
				//case where we have to rotate counterclockwise
				alpha = atan((y1-y2)/(x2-x1));
	//			printf("Rotate %f degrees counterclockwise \n",alpha/(2*PI)*360);
				m1 = cos(alpha);
				m2 = -sin(alpha);
				m3 = sin (alpha);
				m4 = cos(alpha);
			}
		}
		else if (x1 > x2)
		{
			//case where source.x is bigger then target.x, means that associated face lies below edge initially
			if (y1 <= y2)
			{
				//we have to rotate clockwise, but by difference of angle to 180°
				alpha = PI - atan((y2-y1)/(x1-x2));
			//	printf("Rotate %f degrees clockwise \n",alpha/(2*PI)*360);
				m1 = cos(alpha);
				m2 = sin(alpha);
				m3 = -sin (alpha);
				m4 = cos(alpha);

			}
			else
			{
				//we have to rotate counterclockwise, but by additional angle of 180°
				alpha = PI - atan((y1-y2)/(x1-x2));
	//			printf( "Rotate %f degrees counterclockwise \n",alpha/(2*PI)*360);
				m1 = cos(alpha);
				m2 = -sin(alpha);
				m3 = sin (alpha);
				m4 = cos(alpha);
			}
		}

		//apply rotation to first edge
		xx1 = x1*m1 + y1*m2;
		yy1 = x1*m3 + y1*m4;
		xx2 = x2*m1 + y2*m2;
		yy2 = x2*m3 + y2*m4;

		//compute center of rotated hhan
		xc = (xx2-xx1)/2 + xx1;
		yc = (yy2-yy1)/2 + yy1;

		//apply translation such that pC is origin
		xxx1 = xx1 - xc;
		yyy1 = yy1 - yc;
		xxx2 = xx2 - xc;
		yyy2 = yy2 - yc;

		h = 10000000;

		//case where the face into which we bend the edge is the unbounded face (no skeleton region exists)
		Conic_halfedge_handle hithand;
		if (fhan->is_unbounded())
		{
			if (debug)
			{
				printf(" unbounded face \n");
			}

			if (skeletons.find((*indPrimalFaces)[fhan])== skeletons.end())
			{
				if (debug)
				{
					printf("skeleton doesn't exist \n");
				}

				// make polygon
				My_polygon_K outer = getBoundingPolygon();
				My_polygon_with_holes_K poly(outer);
				My_polygon_K hole;


				Halfedge_handle* holes;
				Halfedge_handle hhan;
				for (int i = 0; i < numberOfContinents; i++)
				{
					poly.add_hole(getContinentPolygon(i));

				}
				straightSkeleton = CGAL::create_interior_straight_skeleton_2(poly);
				skeletons[(*indPrimalFaces)[fhan]] = straightSkeleton;
			}
			else
			{
				if (debug)
				{
					printf("skeleton does exist \n");
				}
				straightSkeleton = skeletons[(*indPrimalFaces)[fhan]];
			}

		}

		//case where we bend into an inner face
		else
		{
			if (debug)
			{
				printf(" bounded face \n");
			}

			if (skeletons.find((*indPrimalFaces)[fhan])== skeletons.end())
			{
				if (debug)
				{
					printf("skeleton doesnt exist \n");
				}
				pgn = extractPolygonKFromFace(fhan);
				straightSkeleton = CGAL::create_interior_straight_skeleton_2(pgn.vertices_begin(), pgn.vertices_end());

				skeletons[(*indPrimalFaces)[fhan]] = straightSkeleton;
			}
			else
			{
				if (debug)
				{
					printf("skeleton does exist \n");
				}
				straightSkeleton = skeletons[(*indPrimalFaces)[fhan]];
			}

		}

		//iterate over all skeleton edges belonging to the skeleton face of the edge to be bent (this edge must be the contour edge of the skeleton face)
		Halfedge_iterator skelhit;
		for (  skelhit = (*straightSkeleton).halfedges_begin(); (found == false) && (skelhit != (*straightSkeleton).halfedges_end()); ++skelhit )
		{
			//find corresponding skeleton contour edge
			if (!(skelhit->is_bisector()) && compareEdges(hhan,skelhit))
			{
				if (debug)
				{
					printf("found matching bisector \n");
				}
				found = true;
				skeletonFace = skelhit->face();
				//now traverse all skeleton edges that border the face belonging to the skeleton contour edge hhan (other than actual border contour edge)
				hprov = 10000000;
				for ( skeletonFaceIterator = skeletonFace->halfedge()->next(); skeletonFaceIterator != skeletonFace->halfedge(); skeletonFaceIterator = skeletonFaceIterator->next())
				{
					a = 0;
					m = 0;
					n = 0;

					x3 = CGAL::to_double(skeletonFaceIterator->vertex()->point().x());
					y3 = CGAL::to_double(skeletonFaceIterator->vertex()->point().y());
					x4 = CGAL::to_double(skeletonFaceIterator->opposite()->vertex()->point().x());
					y4 = CGAL::to_double(skeletonFaceIterator->opposite()->vertex()->point().y());

					if (debug)
					{
						printf("p1 = (%f,%f) p2 = (%f,%f) \n",x1,y1,x2,y2);
						printf("p3 = (%f,%f) p4 = (%f,%f) \n",x3,y3,x4,y4);
					}

					//apply rotation to second edge
					xx3 = x3*m1 + y3*m2;
					yy3 = x3*m3 + y3*m4;
					xx4 = x4*m1 + y4*m2;
					yy4 = x4*m3 + y4*m4;

					if (debug)
					{
						printf("p1' = (%f,%f) p2' = (%f,%f) \n",xx1,yy1,xx2,yy2);
						printf("p3' = (%f,%f) p4' = (%f,%f) \n",xx3,yy3,xx4,yy4);
					}

					//apply translation to second edge
					xxx3 = xx3 - xc;
					yyy3 = yy3 - yc;
					xxx4 = xx4 - xc;
					yyy4 = yy4 - yc;


					//let second edge go from left to right too
					if (xxx3 > xxx4)
					{
						//swap p3'' and p4''

						aux1 = xxx3;
						aux2 = yyy3;
						xxx3 = xxx4;
						yyy3 = yyy4;
						xxx4 = aux1;
						yyy4 = aux2;
					}

					if (debug)
					{
						printf("p1'' = (%f,%f) p2'' = (%f,%f) \n",xxx1,yyy1,xxx2,yyy2);
						printf("p3'' = (%f,%f) p4'' = (%f,%f) \n",xxx3,yyy3,xxx4,yyy4);
					}

					a = xxx2;

					//now we are in a simplified setting where the first segment lies on the x-axis and goes from x = -a to x = a
					//the center of the circle lies thus on the y-axis with y = b
					//we want to bend the edge upwards
					//we distinguish several cases wrt the position of the second segment


					if (yyy3 <= 0 && yyy4 <= 0)
					{
						//case1: the circular segment and the line will never intersect, so we set the height to a sufficiently high value
						if (debug)
						{
							printf("case 1 \n");
						}
						hprov = 10000000;
					}
					else if (skelhit->vertex()->point() == skeletonFaceIterator->vertex()->point() )
					{
						//case 0: first and second segment have one endpoint in common
						if (debug)
						{
							printf("case 0.1 \n");
						}
						//calculate the angle between the two segments

						beta = calculateAngle(skelhit,skeletonFaceIterator,skelhit->vertex()->point());
						if (debug)
						{
							printf("angle of %f \n",beta/(2*PI)*360);
						}

						//apply formula for circular segment height (base angle theta = 2*alpha)
						hprov = s/2*tan(beta/2);
					}
					else if (skelhit->opposite()->vertex()->point() == skeletonFaceIterator->vertex()->point() )
					{
						//case 0: first and second segment have one endpoint in common
						if (debug)
						{
							printf("case 0.2 \n");
						}
						//calculate the angle between the two segments

						beta = calculateAngle(skelhit,skeletonFaceIterator,skelhit->opposite()->vertex()->point());
						if (debug)
						{
							printf("angle of %f \n",beta/(2*PI)*360);
						}

						//apply formula for circular segment height (base angle theta = 2*alpha)
						hprov = s/2*tan(beta/2);
					}
					else if (skelhit->opposite()->vertex()->point() == skeletonFaceIterator->opposite()->vertex()->point() )
					{
						//case 0: first and second segment have one endpoint in common
					//	printf("case 0.3 \n");
						//calculate the angle between the two segments

						beta = calculateAngle(skelhit,skeletonFaceIterator,skelhit->opposite()->vertex()->point());
						if (debug)
						{
							printf("angle of %f \n",beta/(2*PI)*360);
						}

						//apply formula for circular segment height (base angle theta = 2*alpha)
						hprov = s/2*tan(beta/2);
					}
					else if (skelhit->vertex()->point() == skeletonFaceIterator->opposite()->vertex()->point() )
					{
						//case 0: first and second segment have one endpoint in common
						if (debug)
						{
							printf("case 0.4 \n");
						}
						//calculate the angle between the two segments

						beta = calculateAngle(skelhit,skeletonFaceIterator,skelhit->vertex()->point());
						if (debug)
						{
							printf("angle of %f \n",beta/(2*PI)*360);
						}

						//apply formula for circular segment height (base angle theta = 2*alpha)
						hprov = s/2*tan(beta/2);
					}
					else if (xxx3 == xxx4)
					{
						//perpendicular cases
						if (xxx3 < xxx2 && xxx3 > xxx1)
						{
							//case2: second edge lies above hhan
							//work on lower vertex of second edge
							if (debug)
							{
								printf("case 2 \n");
							}
							if (yyy3 <= yyy4)
							{
								b1 = (pow(xxx3,2)+pow(yyy3,2)-pow(a,2))/(2*yyy3);
								hprov = sqrt(pow(b1,2)+pow(a,2)) + b1;
							}
							else
							{
								b1 = (pow(xxx4,2)+pow(yyy4,2)-pow(a,2))/(2*yyy4);
								hprov = sqrt(pow(b1,2)+pow(a,2)) + b1;
							}
						}
						else
						{
							//case3: second egde lies to the right or left of hhan
							if (debug)
							{
								printf("case 3 \n");
							}

							b1 = +sqrt(pow(xxx3,2) - pow(xxx2,2));
							b2 = -sqrt(pow(xxx3,2) - pow(xxx2,2)); //can be discarded since intersection lies below hhan (not in circular segment)
							if (debug)
							{
								printf("b1=%f  b2=%f \n",b1,b2);
							}
							if (yyy3 <= yyy4)
							{
								if (b1 >= yyy3 && b1 <= yyy4)
								{
									//cuts segment
									if (debug)
									{
										printf("cuts segment \n");
									}
									hprov = sqrt(pow(b1,2)+pow(a,2)) + b1;
								}
								else
								{
									//cuts line
									if (debug)
									{
										printf("cuts line \n");
									}

									if (yyy3 > 0)
									{
										b1 = (pow(xxx3,2)+pow(yyy3,2)-pow(a,2))/(2*yyy3);
										hprov = sqrt(pow(b1,2)+pow(a,2)) + b1;
									}
									if (yyy4 >0)
									{
										b2 = (pow(xxx4,2)+pow(yyy4,2)-pow(a,2))/(2*yyy4);
										hprov = min(sqrt(pow(b2,2)+pow(a,2)) + b2,hprov);
									}
									if (debug)
									{
										printf("b1=%f  b2=%f \n",b1,b2);
									}
								}
							}
							else
							{
								if (b1 >= yyy4 && b1 <= yyy3)
								{
									//cuts segment
									if (debug)
									{
										printf("cuts segment \n");
									}
									hprov = sqrt(pow(b1,2)+pow(a,2)) + b1;
								}
								else
								{
									//cuts line
									if (debug)
									{
										printf("cuts line \n");
									}

									if (yyy3 > 0)
									{
										b1 = (pow(xxx3,2)+pow(yyy3,2)-pow(a,2))/(2*yyy3);
										hprov = sqrt(pow(b1,2)+pow(a,2)) + b1;
									}
									if (yyy4 >0)
									{
										b2 = (pow(xxx4,2)+pow(yyy4,2)-pow(a,2))/(2*yyy4);
										hprov = min(sqrt(pow(b2,2)+pow(a,2)) + b2,hprov);
									}
									if (debug)
									{
										printf("b1=%f  b2=%f \n",b1,b2);
									}
								}
							}
						}
					}
					else
					{
						//non perpendicular cases (4 and 5)
						m = (yyy4-yyy3)/(xxx4-xxx3);
						n = yyy3 - xxx3*m;

						if (debug)
						{
							printf ("a=%f m=%f n=%f \n", a,m,n);
						}

						if (m==0)
						{
							//parallel case
							if (debug)
							{
								printf("parallel case 5 \n");
							}
							b1 = (pow(n,2)-pow(a,2))/(2*n);
							b2 = -1000000;
							xres1 = 0;
							xres2 = 0;
						}
						else
						{
							if (debug)
							{
								printf("non-parallel case 4 \n");
							}

							b1 = (-n + sqrt(abs(pow(n,2)) - abs(pow(m,4)*pow(a,2)) + abs(pow(m,2)*pow(n,2)) - abs(pow(m,2)*pow(a,2))))/pow(m,2);
							b2 = (-n - sqrt(abs(pow(n,2)) - abs(pow(m,4)*pow(a,2)) + abs(pow(m,2)*pow(n,2)) - abs(pow(m,2)*pow(a,2))))/pow(m,2);
							xres1 = - (m*n - b1*m)/(1+pow(m,2));
							xres2 = - (m*n - b2*m)/(1+pow(m,2));
						}

						if (debug)
						{
							printf("b1=%f  b2=%f \n",b1,b2);
							printf("x1=%f  x2=%f \n",xres1,xres2);
						}


						//one of the solutions intersects above and one below hhan
						if (xres1*m+n >= 0)
						{
							//look if circle really cuts segment or just the line
							if (xres1 >= xxx3 && xres1 <= xxx4)
							{
								//cuts segment
								if (debug)
								{
									printf("cuts segment \n");
								}
								hprov = sqrt(pow(b1,2)+pow(a,2)) + b1;
							}
							else
							{
								//cuts line
								if (debug)
								{
									printf("cuts line \n");
								}

								if (yyy3 > 0)
								{
									b1 = (pow(xxx3,2)+pow(yyy3,2)-pow(a,2))/(2*yyy3);
									hprov = sqrt(pow(b1,2)+pow(a,2)) + b1;
								}
								if (yyy4 >0)
								{
									b2 = (pow(xxx4,2)+pow(yyy4,2)-pow(a,2))/(2*yyy4);
									hprov = min(sqrt(pow(b2,2)+pow(a,2)) + b2,hprov);
								}
							}
						}
						else
						{
							//look if circle really cuts segment or just the line
							if (xres2 >= xxx3 && xres2 <= xxx4)
							{
								//cuts segment
								if (debug)
								{
									printf("cuts segment \n");
								}

								hprov = sqrt(pow(b2,2)+pow(a,2)) + b2;

							}
							else
							{
								//cuts line
								if (debug)
								{
									printf("cuts line \n");
								}

								if (yyy3 > 0)
								{
									b1 = (pow(xxx3,2)+pow(yyy3,2)-pow(a,2))/(2*yyy3);
									hprov = sqrt(pow(b1,2)+pow(a,2)) + b1;
								}
								if (yyy4 >0)
								{
									b2 = (pow(xxx4,2)+pow(yyy4,2)-pow(a,2))/(2*yyy4);
									hprov = min(sqrt(pow(b2,2)+pow(a,2)) + b2,hprov);
								}
							}
						}
					}


					if (debug)
					{
						printf("hprov=%f \n",hprov);
					}
					if (hprov < h )
					{
						h = hprov;
					}
				}
			}
		}

		if (debug)
		{
			printf ("got height of %f \n", h);
		}

		A = (0.5 * atan(2*h/s) * pow(4*pow(h,2) + pow(s,2),2) + h*s*(4*pow(h,2) - pow(s,2))) / (16*pow(h,2));


		return A;

}


/*
 * computes if the target vertex of hhan (at the outer face) is reflex where p is the other vertex on the adjacent egde
 */
bool MyWindow::vertexIsReflex(Conic_halfedge_handle hhan, Arr_conic_point_2 p)
{

	double x1 = CGAL::to_double(hhan->source()->point().x());
	double y1 = CGAL::to_double(hhan->source()->point().y());
	double x2 = CGAL::to_double(hhan->target()->point().x());
	double y2 = CGAL::to_double(hhan->target()->point().y());
	double x3 = CGAL::to_double(p.x());
	double y3 = CGAL::to_double(p.y());
	double s =  sqrt(pow(x1 - x2,2)  + pow(y1-y2,2));
	double m1,m2,m3,m4,alpha,xx1,xx2,xx3,yy1,yy2,yy3;


	// get rotation matrix such that after its application the edge is horizontal, points from left to right and the associated face lies above edge (to its left)
	if (x1 == x2)
	{
		//vertical case
		if (y1 > y2)
		{
			//case where we have to rotate clockwise
			alpha = PI/2;;
	//		printf("Rotate %f degrees counterclockwise \n",alpha/(2*PI)*360);
			m1 = cos(alpha);
			m2 = -sin(alpha);
			m3 = sin (alpha);
			m4 = cos(alpha);
		}
		else
		{
			//case where we have to rotate counterclockwise
			alpha = PI/2;
	//		printf("Rotate %f degrees clockwise \n",alpha/(2*PI)*360);
			m1 = cos(alpha);
			m2 = sin(alpha);
			m3 = -sin (alpha);
			m4 = cos(alpha);
		}
	}
	else if (x1 < x2)
	{
		//case where source.x is smaller then target.x, means that associated face lies above edge initially
		if (y1 <= y2)
		{
			//case where we have to rotate clockwise
			alpha = atan((y2-y1)/(x2-x1));
	//		printf("Rotate %f degrees clockwise \n",alpha/(2*PI)*360);
			m1 = cos(alpha);
			m2 = sin(alpha);
			m3 = -sin (alpha);
			m4 = cos(alpha);
		}
		else
		{
			//case where we have to rotate counterclockwise
			alpha = atan((y1-y2)/(x2-x1));
	//		printf("Rotate %f degrees counterclockwise \n",alpha/(2*PI)*360);
			m1 = cos(alpha);
			m2 = -sin(alpha);
			m3 = sin (alpha);
			m4 = cos(alpha);
		}
	}
	else if (x1 > x2)
	{
		//case where source.x is bigger then target.x, means that associated face lies below edge initially
		if (y1 <= y2)
		{
			//we have to rotate clockwise, but by difference of angle to 180°
			alpha = PI - atan((y2-y1)/(x1-x2));
		//	printf("Rotate %f degrees clockwise \n",alpha/(2*PI)*360);
			m1 = cos(alpha);
			m2 = sin(alpha);
			m3 = -sin (alpha);
			m4 = cos(alpha);

		}
		else
		{
			//we have to rotate counterclockwise, but by additional angle of 180°
			alpha = PI - atan((y1-y2)/(x1-x2));
	//		printf( "Rotate %f degrees counterclockwise \n",alpha/(2*PI)*360);
			m1 = cos(alpha);
			m2 = -sin(alpha);
			m3 = sin (alpha);
			m4 = cos(alpha);
		}
	}

	//apply rotation to hhan, afterwards hhan lies horizontally and points from left to right (outer face is above hhan)
	xx1 = x1*m1 + y1*m2;
	yy1 = x1*m3 + y1*m4;
	xx2 = x2*m1 + y2*m2;
	yy2 = x2*m3 + y2*m4;

	//apply rotation to p
	xx3 = x3*m1 + y3*m2;
	yy3 = x3*m3 + y3*m4;

	//if the y value of rotated p is bigger than the y value of horizontally rotated hhan, then the vertex is reflex
	return (yy3 > yy2);
}


//get the angle at shared vertex between the two edges corresponding to the halfedges hh1 and hh2
//precondition: hh1 and hh2 share point p1
double MyWindow::calculateAngle(Halfedge_handle hh1, Halfedge_handle hh2, K_point p1)
{

	//get the other two vertices (one vertex of hh1, one of hh2, those that don't coincide with p1)
	K_point p2;
	K_point p3;
	if (p1 != hh1->vertex()->point())
	{
		 p2 = hh1->vertex()->point();
	}
	else
	{
		 p2 = hh1->opposite()->vertex()->point();
	}
	if (p1 != hh2->vertex()->point())
	{
		 p3 = hh2->vertex()->point();
	}
	else
	{
		 p3 = hh2->opposite()->vertex()->point();
	}

	//get the lengths of all three triangle sides
	double b=sqrt(fabs(p1.x() - p2.x())*fabs(p1.x() - p2.x())+
			fabs(p1.y() - p2.y())*fabs(p1.y() - p2.y()));
	double c=sqrt(fabs(p1.x() - p3.x())*fabs(p1.x() - p3.x())+
				fabs(p1.y() - p3.y())*fabs(p1.y() - p3.y()));
	double a=sqrt(fabs(p2.x() - p3.x())*fabs(p2.x() - p3.x())+
				fabs(p2.y() - p3.y())*fabs(p2.y() - p3.y()));

	//apply angle formula
	double alpha = acos((b*b+c*c-a*a)/(2*b*c));

	return alpha;
}


//get the angle at shared vertex between the two edges corresponding to the halfedges hh1 and hh2
//precondition: hh1 and hh2 share point p1
double MyWindow::calculateAngle(Conic_halfedge_handle hh1, Conic_halfedge_handle hh2, Arr_conic_point_2 p1)
{

	//get the other two vertices (one vertex of hh1, one of hh2, those that don't coincide with p1)
	Arr_conic_point_2 p2;
	Arr_conic_point_2 p3;

	if (p1 != hh1->source()->point())
	{
		 p2 = hh1->source()->point();
	}
	else
	{
		 p2 = hh1->target()->point();
	}
	if (p1 != hh2->source()->point())
	{
		 p3 = hh2->source()->point();
	}
	else
	{
		 p3 = hh2->target()->point();
	}

	//get the lengths of all three triangle sides
	double b=sqrt(pow(CGAL::to_double(p1.x()) - CGAL::to_double(p2.x()),2) +
				  pow(CGAL::to_double(p1.y()) - CGAL::to_double(p2.y()),2));
	double c=sqrt(pow(CGAL::to_double(p1.x()) - CGAL::to_double(p3.x()),2) +
				  pow(CGAL::to_double(p1.y()) - CGAL::to_double(p3.y()),2));
	double a=sqrt(pow(CGAL::to_double(p2.x()) - CGAL::to_double(p3.x()),2) +
				  pow(CGAL::to_double(p2.y()) - CGAL::to_double(p3.y()),2));

	//apply angle formula
	double alpha = acos((b*b+c*c-a*a)/(2*b*c));

	return alpha;
}



//checks whether the two egdes of different types are identical (have same coordinates
//optionally it can be abstracted from the direction of the edges (commented part at the end of return statement)
bool MyWindow::compareEdges(Conic_halfedge_handle  h1,Halfedge_handle  h2)
{
	K_point source (CGAL::to_double(h1->source()->point().x()),CGAL::to_double(h1->source()->point().y()));
	K_point target (CGAL::to_double(h1->target()->point().x()),CGAL::to_double(h1->target()->point().y()));
	//check for equality of the two half-edges
	//maybe need to consider possible different orientations of edges?? if so, uncomment last two lignes

	return  source == h2->opposite()->vertex()->point() &&
			target == h2->vertex()->point();
			/*||
			(source == h2->vertex()->point() &&
			target == h2->opposite()->vertex()->point()));
			*/
}




//calculates the height of a circular arc on straight line segment of length s where the area of the circular segment is parameter a
//for the moment this uses only an approximation (binary search of suitable height with 1% error)! exact computation has to be implemented (solving circular segment equation)
double MyWindow::areaToArcHeight(double a, Arr_conic_point_2 p1, Arr_conic_point_2 p2)
{
	//get segment length
	double s = sqrt(pow(CGAL::to_double(p1.x()) - CGAL::to_double(p2.x()),2) + pow(CGAL::to_double(p1.y()) - CGAL::to_double(p2.y()),2));


	//start binary search with a height that implies a half-circle
	double h = 0.1;
	double area = (0.5*atan(2*h/s)*(4*h*h+s*s)*(4*h*h+s*s)+h*s*(4*h*h-s*s)) / (16*h*h);

	while (area < a)
	{
		h = h*1.01;
		area = (0.5*atan(2*h/s)*(4*h*h+s*s)*(4*h*h+s*s)+h*s*(4*h*h-s*s)) / (16*h*h);
	}

	return h/1.01;
}
