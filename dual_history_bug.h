#ifndef DUAL_HISTORY_BUG
#define DUAL_HISTORY_BUG

namespace CGAL {


template <class GeomTraits_, class TopTraits_>
class Dual<Arrangement_on_surface_with_history_2<GeomTraits_,TopTraits_> >
{
public:
  
  typedef GeomTraits_                          Geometry_traits_2;
  typedef TopTraits_                           Topology_traits;
  typedef CGAL::Arrangement_on_surface_with_history_2<Geometry_traits_2, Topology_traits>
                                               Arrangement_on_surface_with_history_2;

  typedef typename Arrangement_on_surface_with_history_2::Size              Size;
  typedef typename Arrangement_on_surface_with_history_2::Face_handle       Vertex_handle;
  typedef typename Arrangement_on_surface_with_history_2::Halfedge_handle   Edge_handle;

  typedef typename Arrangement_on_surface_with_history_2::Face_iterator     Vertex_iterator;
  typedef typename Arrangement_on_surface_with_history_2::Halfedge_iterator Edge_iterator;

protected:

  typedef typename Arrangement_on_surface_with_history_2::Face_handle       Face_handle;
  typedef typename Arrangement_on_surface_with_history_2::Ccb_halfedge_circulator
                                                      Ccb_halfedge_circulator;
  typedef typename Arrangement_on_surface_with_history_2::Outer_ccb_iterator
                                                      Outer_ccb_iterator;
  typedef typename Arrangement_on_surface_with_history_2::Inner_ccb_iterator
                                                      Inner_ccb_iterator;

  /*! \class
   * Iterator over the neighbors of a dual vertex (a face in the primal
   * arrangement).
   * These neighbors are the adjacent faces along the outer boundaries of the
   * face and its inner boundaries.
   */
  class Face_neighbor_iterator
  {
    typedef Face_neighbor_iterator               Self;

  public:

    typedef std::forward_iterator_tag            iterator_category;
    typedef Edge_handle                          value_type;
    typedef value_type                           reference;
    typedef value_type*                          pointer;
    typedef int                                  difference_type;

  private:

    Outer_ccb_iterator       _outer_ccb_iter;
    Inner_ccb_iterator       _inner_ccb_iter;
    Ccb_halfedge_circulator  _ccb_curr;
    Ccb_halfedge_circulator  _ccb_first;
    Face_handle              _face;
    bool                     _out;
    Edge_handle              _hh;
    bool                     _end;

  public:

    /*! Default constructor. */
    Face_neighbor_iterator () :
      _end (true)
    {}

    /*!
     * Constructor.
     * \param face The face (dual vertex).
     * \param out_edges Do we need the outgoing or the ingoing halfedges.
     * \param start Should we start traversing the edges.
     *              If false, we construct a past-the-end iterator.
     */
    Face_neighbor_iterator (Face_handle face, 
                            bool out_edges,
                            bool start) :
      _face (face),
      _out (out_edges),
      _end (! start)
    {
      CGAL_precondition (! face->is_fictitious());

      if (start)
      {
        _outer_ccb_iter = _face->outer_ccbs_begin();
        _inner_ccb_iter = _face->inner_ccbs_begin();

        if (_outer_ccb_iter != _face->outer_ccbs_end())
        {
          // Start from the first outer CCB, if one exists.
          _ccb_curr = _ccb_first = *_outer_ccb_iter;
        }
        else if (_inner_ccb_iter != face->inner_ccbs_end())
        {
          // Otherwise, start from the first inner CCB.
          _ccb_curr = _ccb_first = *_inner_ccb_iter;
        }
        else
        {
          // In this case there are no CCBs to traverse:
          _end = true;
          return;
        }

        _hh = this->_dereference();

        // In case the incident face of the twin halfedge is fictitious,
        // skip it and proceed to the next edge.
        if (_hh->is_fictitious())
          ++(*this);
      }
      else // end iterator.
      {
        _outer_ccb_iter = _face->outer_ccbs_end();
        _inner_ccb_iter = _face->inner_ccbs_end();
      }
    }  

    /*! Equality operators. */
    bool operator== (const Self& it) const
    {
      return (this->_equal(it));
    }
    
    bool operator!= (const Self& it) const
    {
      return (! this->_equal(it));
    }
    
    /*! Dereference operators. */
    reference operator* () const
    {
      return (_hh);
    }

    pointer operator-> () const
    {
      return (&_hh);
    }
    
    /* Increment operators. */
    Self& operator++ ()
    {
      do
      {
        this->_increment();
        if (_end)
          return (*this);

        _hh = this->_dereference();

      } while (_hh->is_fictitious());

      return (*this);
    }

    Self operator++ (int )
    {
      Self tmp = *this;

      do
      {
        this->_increment();
        if (_end)
          return (tmp);

        _hh = this->_dereference();

      } while (_hh->is_fictitious());

      return (tmp);
    }

  private:

    /*! Check two iterators for equality. */
    bool _equal (const Self& it) const
    {
      return (_out == it._out && _face == it._face &&
              ((_end && it._end) ||
               (_outer_ccb_iter == it._outer_ccb_iter &&
                _inner_ccb_iter == it._inner_ccb_iter &&
                _ccb_curr == it._ccb_curr)));
    }

    /*! Derefernce the current circulator. */
    Edge_handle _dereference () const
    {
      if (_out)
        return (_ccb_curr);
      else
        return (_ccb_curr->twin());
    }

    // Increments of the iterator.
    void _increment ()
    {
      CGAL_assertion (! _end);

      // If we have not traversed the entire CCB in full, move to the next
      // halfedge along the current CCB.
      ++_ccb_curr;

      if (_ccb_curr != _ccb_first)
        return;

      // In this case we have completed the current CCB and we have to move
      // to the next one.
      if (_outer_ccb_iter != _face->outer_ccbs_end())
      {
        // Try to move to the next outer CCB.
        ++_outer_ccb_iter;
        if (_outer_ccb_iter != _face->outer_ccbs_end())
        {
          _ccb_curr = _ccb_first = *_outer_ccb_iter;
          return;
        }

        // In this case we start traversing the inner CCBs.
        if (_inner_ccb_iter != _face->inner_ccbs_end())
        {
          CGAL_assertion (_inner_ccb_iter == _face->inner_ccbs_begin());

          // Otherwise, start from the first inner CCB.
          _ccb_curr = _ccb_first = *_inner_ccb_iter;
          return;
        }
      }
      else if (_inner_ccb_iter != _face->inner_ccbs_end())
      {

        // In this case we have already traversed all outer CCBs (and at least
        // one inner CCB), so we try to move to the next inner CCB.
        ++_inner_ccb_iter;
        if (_inner_ccb_iter != _face->inner_ccbs_end())
        {
          // Otherwise, start from the first inner CCB.
          _ccb_curr = _ccb_first = *_inner_ccb_iter;
          return;
        }
      }

      // In this case we finished traversing all outer and inner CCBs:
      _end = true;
      return;
    }

  };

  // Data members:
  mutable Arrangement_on_surface_with_history_2    *p_arr;    // The primal arrangement.
  
public:

  typedef Face_neighbor_iterator            Incident_edge_iterator;

  /*! Default constructor. */
  Dual () :
    p_arr (NULL)
  {}

  /*! Constructor from an arrangement. */
  Dual (const Arrangement_on_surface_with_history_2& arr) :
    p_arr (const_cast<Arrangement_on_surface_with_history_2 *> (&arr))
  {}

  /*! Get the primal arrangement (const version). */
  const Arrangement_on_surface_with_history_2* arrangement () const
  {
    return (p_arr);
  }

  /*! Get the primal arrangement (non-const version). */
  Arrangement_on_surface_with_history_2* arrangement ()
  {
    return (p_arr);
  }

  /*! Get the number of vertices (face of the primal arrangement). */
  Size number_of_vertices () const
  {
    return (p_arr->number_of_faces());
  }

  /*! Traverse the vertices (faces of the primal arrangement). */
  Vertex_iterator vertices_begin () const
  {
    return (p_arr->faces_begin());
  }

  Vertex_iterator vertices_end () const
  {
    return (p_arr->faces_end());
  }

  /*! Get the number of edges. */
  Size number_of_edges () const
  {
    return (p_arr->number_of_halfedges());
  }

  /*! Traverse the edges. */
  Edge_iterator edges_begin () const
  {
    return (p_arr->halfedges_begin());
  }

  Edge_iterator edges_end () const
  {
    return (p_arr->halfedges_end());
  }

  /*!
   * Get the dual vertex-degree (number of edges forming the face boundary).
   */
  Size degree (Vertex_handle v) const
  {
    Incident_edge_iterator   begin = Incident_edge_iterator (v, true, true);
    Incident_edge_iterator   end = Incident_edge_iterator (v, false, true);
    Size                     deg = 0;

    while (begin != end)
    {
      deg++;
      ++begin;
    }

    return (deg);
  }

  /*! Traverse the outgoing edges of a given vertex. */
  Incident_edge_iterator out_edges_begin (Vertex_handle v) const
  {
    return (Incident_edge_iterator (v, true, true));
  }

  Incident_edge_iterator out_edges_end (Vertex_handle v) const
  {
    return (Incident_edge_iterator (v, true, false));
  }

  /*! Traverse the ingoing edges of a given vertex. */
  Incident_edge_iterator in_edges_begin (Vertex_handle v) const
  {
    return (Incident_edge_iterator (v, false, true));
  }

  Incident_edge_iterator in_edges_end (Vertex_handle v) const
  {
    return (Incident_edge_iterator (v, false, false));
  }
};


template <class Traits_, class Dcel_>
class Dual<Arrangement_with_history_2<Traits_, Dcel_> > :
  public Dual<CGAL::Arrangement_on_surface_with_history_2<
             typename CGAL::Arrangement_with_history_2<Traits_, Dcel_>::Geometry_traits_2,
             typename CGAL::Default_planar_topology<Traits_, Dcel_>::Traits> >
{
  typedef Traits_                                                     Traits_2;
  typedef Dcel_                                                       Dcel;


  typedef typename Default_planar_topology<Traits_2, Dcel>::Traits 
  Topology_traits;
  
  typedef Dual<CGAL::Arrangement_on_surface_with_history_2<
                 typename CGAL::Arrangement_with_history_2<Traits_2, Dcel>::Geometry_traits_2,
                 Topology_traits> >  Base;

public:

  /*! Default constructor. */
  Dual () :
    Base()
  {}

  /*! Constructor from an arrangement. */
  Dual (const CGAL::Arrangement_with_history_2<Traits_2, Dcel>& arr) :
      Base (arr)
  {  }
};
}

#endif // DUAL_HISTORY_BUG
