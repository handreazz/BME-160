import numpy
class Triad :
    """
    Calculate angles and distances among a triad of points.
 
    Author: Eric Mockler
    Date: April 14, 2019
    Points can be supplied in any dimensional space as long as they are consistent.
    Points are supplied as tupels in n-dimensions, and there should be three
    of those to make the triad. Each point is positionally named as p,q,r
    and the corresponding angles are then angleP, angleQ and angleR.
    Distances are given by dPQ(), dPR() and dQR()
 
    Required Modules: NumPy
    initialized: 3 positional tuples converted to NumPy arrays, representing Points in n-space
             p1 = Triad( p=(1,0,0), q=(0,0,0), r=(0,1,0) )
    attributes: p,q,r the 3 tuples representing points in N-space
    methods:  angleP(), angleR(), angleQ() angles measured in radians
          dPQ(), dPR(), dQR() distances in the same units of p,q,r
 
    """
 
    def __init__(self,p,q,r) :
        """ Construct a Triad. 
        
        Example construction:
            p1 = Triad( p=(1.,0.,0.), q=(0.,0.,0.), r=(0.,0.,0.) ). 
        """
        self.p = numpy.asarray(p)
        self.q = numpy.asarray(q)
        self.r = numpy.asarray(r)
    
# calculate lengths(distances) of segments PQ, PR and QR
    def dPQ (self):
        """ Provides the distance between point p and point q """
        return numpy.linalg.norm(self.p-self.q)
    def dPR (self):
        """ Provides the distance between point p and point r """
        return numpy.linalg.norm(self.p-self.r)
    def dQR (self):
        """ Provides the distance between point q and point r """
        return numpy.linalg.norm(self.q-self.r)

    def angleP (self) :
        """ Provides the angle made at point p by segments pq and pr (radians). """
        return numpy.arccos(numpy.dot((self.q-self.p), (self.r-self.p))/(numpy.linalg.norm(self.q-self.p) * numpy.linalg.norm(self.r-self.p)))
    def angleQ (self) :
        """ Provides the angle made at point q by segments qp and qr (radians). """
        return numpy.arccos(numpy.dot((self.p-self.q), (self.r-self.q))/(numpy.linalg.norm(self.p-self.q) * numpy.linalg.norm(self.r-self.q)))
    def angleR (self) :
        """ Provides the angle made at point r by segments rp and rq (radians). """
        return numpy.arccos(numpy.dot((self.p-self.r), (self.q-self.r))/(numpy.linalg.norm(self.p-self.r) * numpy.linalg.norm(self.q-self.r)))
