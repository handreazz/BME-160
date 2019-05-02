#!/usr/bin/env python3 
# Name: Eric Mockler (emockler) 

'''
Calculates C-N-Ca bond angles and legnths given the atoms' 3D coordinates in space

in:
    C = (39.447, 94.657, 11.824) N = (39.292, 95.716, 11.027) Ca = (39.462, 97.101, 11.465) 
out:
    N-C bond length = 1.33
    N-Ca bond length = 1.46
    C-N-Ca bond angle = 124.0
'''

import math
class Triad :
    """
    Calculate angles and distances among a triad of points.
 
    Author: David Bernick
    Date: March 21, 2013
    Points can be supplied in any dimensional space as long as they are consistent.
    Points are supplied as tupels in n-dimensions, and there should be three
    of those to make the triad. Each point is positionally named as p,q,r
    and the corresponding angles are then angleP, angleQ and angleR.
    Distances are given by dPQ(), dPR() and dQR()
 
    Required Modules: math
    initialized: 3 positional tuples representing Points in n-space
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
        self.p = p
        self.q = q
        self.r = r
# private helper methods
    def d2 (self,a,b) : # calculate squared distance of point a to b
        return float(sum((ia-ib)*(ia-ib)  for  ia,ib in zip (a,b)))
    
    def dot (self,a,b) : # dotProd of standard vectors a,b
        return float(sum(ia*ib for ia,ib in zip(a,b)))
    
    def ndot (self,a,b,c) : # dotProd of vec. a,c standardized to b
        return float(sum((ia-ib)*(ic-ib) for ia,ib,ic in zip (a,b,c)))
    
# calculate lengths(distances) of segments PQ, PR and QR
    def dPQ (self):
        """ Provides the distance between point p and point q """
        return math.sqrt(self.d2(self.p,self.q))
    
    def dPR (self):
        """ Provides the distance between point p and point r """
        return math.sqrt(self.d2(self.p,self.r))
    
    def dQR (self):
        """ Provides the distance between point q and point r """
        return math.sqrt(self.d2(self.q,self.r))
    
    def angleP (self) :
        """ Provides the angle made at point p by segments pq and pr (radians). """
        return math.acos(self.ndot(self.q,self.p,self.r) /   math.sqrt(self.d2(self.q,self.p)*self.d2(self.r,self.p)))
    
    def angleQ (self) :
        """ Provides the angle made at point q by segments qp and qr (radians). """
        return math.acos(self.ndot(self.p,self.q,self.r) /  math.sqrt(self.d2(self.p,self.q)*self.d2(self.r,self.q)))
 
    def angleR (self) :
        """ Provides the angle made at point r by segments rp and rq (radians). """
        return math.acos(self.ndot(self.p,self.r,self.q) /  math.sqrt(self.d2(self.p,self.r)*self.d2(self.q,self.r)))

class FindCoords(str) :
    """
    Finds bond lengths and angles from user-inputted C, N, and Ca coordinates
    """
    def _length_(self):
        """ Provides length of input """
        return len(self)
    def _tuple_(self, xTuple):
        """ Returns a tuple containing float items """
        for i, coords in enumerate(xTuple):
            xTuple[i] = float(coords)
        return tuple(xTuple)    

    # end of private helper methods
    def parse(self):
        """ Returns user input into a 3d-coordinate tuple in the order C, N, Ca"""
        inCoords = self
        inCoords = inCoords.upper()
        # remove input that's not a coordinate or useful delimiter
        for c in ('N', 'CA', 'C', '=', ' '):
            inCoords = inCoords.replace(c, '')
        # remove opening parenthesis & use closing parenthesis as delimiter 
        inCoords = inCoords.replace('(', '')
        inCoords = inCoords.split(')')
        # create C 3-tuple
        cTuple = inCoords[0].split(',')
        cTuple = self._tuple_(cTuple)
        # create N 3-tuple
        nTuple = inCoords[1].split(',')
        nTuple = self._tuple_(nTuple)
        # create Ca 3-tuple
        caTuple = inCoords[2].split(',')
        caTuple = self._tuple_(caTuple)
        return tuple([cTuple, nTuple, caTuple])
    
    def printCalc(self, rawCoords):
        """ Given a triad with 3D coordinates of (C, N, Ca), N-C and N-Ca bond lengths & C-N-Ca bond angle are printed"""
        triadPoints = Triad(p = (rawCoords[0]), q = (rawCoords[1]), r = (rawCoords[2]))
        nToC = triadPoints.dPQ()
        nToCa = triadPoints.dQR()
        bondAngle = math.degrees(triadPoints.angleQ())
        print(
"""
N-C bond length = {0:.2f}
N-Ca bond length = {1:.2f}
C-N-Ca bond angle = {2:.1f}
""".format(nToC, nToCa, bondAngle))
def main():
    """ Instructs user input of C, N, and Ca coordinates, and prints calculated bond lengths and angles"""
    rawCoords = input(
"""
Enter xyz coordinates in the format:
C = (x1, y1, z1) N = (x2, y2, z2) Ca = (x3, y3, z3)
""" )    
    rawCoords = FindCoords(rawCoords)
    parsedCoords = rawCoords.parse()
    rawCoords.printCalc(parsedCoords)
main()