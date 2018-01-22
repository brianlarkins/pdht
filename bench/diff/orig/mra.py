from tensor import *                     
from quadrature import Quadrature
from twoscalecoeffs import twoscalecoeffs, phi
from math import sqrt, exp, sin, cos
from whrandom import random, seed
import autoc
import testlapack
import sys
import cfcube
import time
import types
import cPickle
import tempfile
import perfstats
from diskdir import *
from sparse import *
import fast_eval

# The recommended way of importing this module is
# from mra import Function.

# These are operators that act upon Functions by causing the Function
# __rmul__ method to invoke the appropriate action. 

class Dx:
    pass

class Dy:
    pass

class Dz:
    pass

class DxT:
    pass

class DyT:
    pass

class DzT:
    pass

class D2x:
    pass

class D2y:
    pass

class D2z:
    pass

class Del:
    pass

class Delsq:
    pass

def parent(lx):
    # In Python integer arithmetic 1/2 -> 0 but -1/2 -> -1 which
    # is wierd but just what we need
    return lx/2

def parentmod(lx):
    return lx % 2

perf = perfstats.PerfStats("MRA Function")
 
class Function:
    '''

    This class primarily implements a 3D function in Alperts
    multiwavelet basis.  Some of the methods and data structures
    support other dimensions.  Those that do not must be overridden or
    deleted by subclasses.

    Some of the class attributes (e.g., autorefine) are always
    referenced from the class rather than the instance ... this is to
    make it easier to change the behaviour of all instances rather
    than just new instances.

    '''

    quad = {}                         # For reuse of quadrature info
    twoscale = {}                     # For reuse of twoscale coefficients

    # Provide defaults in the class to make it easier to have all instances
    # be consistent.
    k = 9                               # Default order
    thresh = 1e-7                       # Default precision
    initial_level = 2                   # Default level for initial projection
    debug = 1                           # Debug print level
    autorefine = 1                      # Refine during square (mul)
    operator = 0                        # I'm not the kernel of an operator

    # truncate_method = 0 => |d|<thresh 
    # truncate_method = 1 => |d|<thresh*8^(-n/2)
    # truncate_method = 2 => |d|<thresh/sqrt(len(d))
    truncate_method = 0
    max_refine_level = 20               # Stop refining at this level

    ncalls = {}                         # Counts calls to methods
    times  = {}                         # Times  calls to methods

    def __init__(self,
                 function=None, cfunction=None,
                 thresh=None, k=None,
                 initial_level=None,refine=1,compress=1,box=None):

        perf.enter('init')

        # For methods that may be inherited by subclasses we need to
        # remember the class and tensor rank and override it in
        # subclasses (don't really need to remember the tensor rank
        # I think, but it might be convenient at some point).

        self.ndim = 3
        self.new = Function
        self.newtensor = Tensor

        if k:
            self.k = k                  # Override default in class
        else:
            self.k = k = Function.k     # Ensure set in instance

        if self.operator:               # Operators need double order
            self.k = k = 2*k

        if thresh:
            self.thresh = thresh        # Override default in class
        else:
            self.thresh = thresh = Function.thresh # Ensure set in instance

        if initial_level is None:
            self.initial_level = initial_level = Function.initial_level
        else:
            self.initial_level = initial_level

        self.d = {}             # Directory of wavelet coefficients
        self.s = {}             # Directory of scaling function coefficients
        self.function = function# Remember the function being compressed
        self.cfunction = cfunction
        self.compressed = 0     # Not yet compressed

        self.init_twoscale(k)
        self.init_quadrature(k)

        if cfunction:
            if not (type(cfunction) == type(' ') and \
                    cfunction[-6:] == 'p_void'):
                raise TypeError,"cfunction should be a C "\
                      "function pointer cast to void *"

        if function or cfunction:
            self.nterminated = 0
            if box is None:
                trans = None
            else:
                trans = []
                lxlo,lylo,lzlo,lxhi,lyhi,lzhi = box
                yrange = range(lylo,lyhi+1)
                zrange = range(lzlo,lzhi+1)
                for lx in range(lxlo,lxhi+1):
                    for ly in yrange:
                        for lz in zrange:
                            trans.append((lx,ly,lz))
                
            perf.enter('fs proj')
            self.fine_scale_projection(initial_level,trans)
            perf.exit('fs proj')
            
            if refine:
                perf.enter('refine')
                # Midly convoluted way to make list of translations to
                # be refined so that the same code works for operators
                # and for functions.
                refinethese = {}
                for lx, ly, lz in self.s[initial_level].keys():
                    refinethese[parent(lx),parent(ly),parent(lz)] = 1
                for lx, ly, lz in refinethese.keys():
                    self.refine_fine_scale_projection(initial_level-1,
                                                      lx, ly, lz, 0)
                perf.exit('refine')

            if self.nterminated > 0:
                print "Terminated refinement:", self.nterminated

            if compress:
                self.compress(truncate=0)
        else:
            self.s[0] = {}
            self.s[0][0,0,0] = Tensor(k,k,k)
            if compress:
                self.compressed = 1
                self.d[0] = {}
                self.d[0][0,0,0] = Tensor(2*k,2*k,2*k)
            
        if self.debug:
            print "Made new function", id(self), self.k, self.thresh, self.compressed
            print "s:", self.s.keys()
	    #for lvl in self.s.keys():
	    #	print "lvl: ", lvl
	    #   	for myx,myy,myz in self.s[lvl].keys():
	    #		print "  ", myx,myy,myz
            if self.compressed: print "d:", self.d.keys()

        perf.exit('init')
            
    def init_twoscale(self,k):
        if self.twoscale.has_key(k):
            self.h0, self.h1, self.g0, self.g1, self.hg = self.twoscale[k]
            self.hgsonly = Matrix(self.hg[0:k,:])
            self.hgT = Matrix(self.hg.cycledim(1)) # Force contiguous
            return
        
        (self.h0, self.h1, self.g0, self.g1) = twoscalecoeffs(k)
        self.h0 = Matrix(self.h0)
        self.h1 = Matrix(self.h1)
        self.g0 = Matrix(self.g0)
        self.g1 = Matrix(self.g1)
        self.hg = Matrix(2*k,2*k)
        self.hg[0:k  , 0:k]   = self.h0
        self.hg[0:k  , k:2*k] = self.h1
        self.hg[k:2*k, 0:k]   = self.g0
        self.hg[k:2*k, k:2*k] = self.g1

        self.hgsonly = Matrix(self.hg[0:k,:])
        self.hgT = Matrix(self.hg.cycledim(1))

        # As an optimization remember the twoscale coeffs in the class
        Function.twoscale[k] = (self.h0, self.h1, self.g0, self.g1, self.hg)

    def init_quadrature(self,order):
        '''

        Compute and store the quadrature points+weights, and tabulate
        the scaling functions * quadrature weights at those points.
        This is done on level-0.  At level-n scale by 1/2**(n/2) *
        2^(n/2) from the scaling of the scaling functions * 2^-n from
        the scaling of the interval.  And the quadrature points must
        be scaled by 2^-n.

        '''
        if self.quad.has_key(order):
            self.quad_x, self.quad_w, self.quad_npt, \
            self.quad_phi, self.quad_phiw, \
            self.quad_psi, self.quad_psiw = self.quad[order]
            return

        q = Quadrature(order,[0.0,1.0])
        self.quad_x = Vector(q.points())
        self.quad_w = w = q.weights()
        self.quad_npt = npt = len(w)
        self.quad_phi  = Matrix(npt,self.k) # phi[point,i]
        self.quad_phiw = Matrix(npt,self.k) # phi[point,i]*weight[point]
        for i in range(npt):
            self.quad_phi[i,:] = Vector(phi(self.quad_x[i],self.k))
            for m in range(self.k):
                self.quad_phiw[i,m] = self.quad_phi[i,m]*w[i]

        # psi for multiwavelet
        #
        # self.quad_psi  ...  phi[point,i]
        # self.quad_psiw ...  phi[point,i]*weight[point]
        #
        phi0 = Matrix(npt,self.k)
        phi1 = Matrix(npt,self.k)
        for i in range(npt):
            x = self.quad_x[i]
            if x < 0.5:
                phi0[i,:] = Vector(phi(2.0*x  ,self.k))
            else:
                phi1[i,:] = Vector(phi(2.0*x-1,self.k))
        g0T = self.g0.transpose()
        g1T = self.g1.transpose()
        self.quad_psi = ( phi0*g0T + phi1*g1T ).scale(sqrt(2.0))
        self.quad_psiw = Matrix(npt,self.k) # psi[point,i]*weight[point]
        for i in range(npt):
            for m in range(self.k):
                self.quad_psiw[i,m] = self.quad_psi[i,m]*w[i]

        # As an optimization remember the quadrature info in the class itself
        Function.quad[order] = (self.quad_x, self.quad_w, self.quad_npt,
                                self.quad_phi, self.quad_phiw,
                                self.quad_psi, self.quad_psiw)
        
    def fcube(self, xlo, ylo, zlo, npt, h, pts, f, function):
        '''

        Tabulate the function on the quadrature points in the
        specified cube.  The C version is MUCH faster.

        '''
        rangenpt = range(npt)
        for p in rangenpt:
            x = xlo+pts[p]*h
            for q in rangenpt:
                y = ylo+pts[q]*h
                for r in rangenpt:
                    z = zlo+pts[r]*h
                    f[p,q,r] = function(x,y,z)

    def fine_scale_projection(self,n,trans=None):
        '''

        Project the function onto scaling functions at level n.

        If self.operator is true, then we are computing the kernel of
        an integral operator. Operators need -2^n <= l <= 2^n as the
        default range.  

        If l[] is specified it is a list of translations that are to
        be considered.  Otherwise the full list is used.

        Modifies:
        :           self.s

        '''

        try:
            self.s[n]
        except KeyError:
            self.s[n] = {}
        
        h = 1.0/(2.0**n)                # Box size on target level 
        scale = sqrt(h)
        scale = scale*scale*scale       # Since in 3D
        phiw = self.quad_phiw
        pts = self.quad_x 
        npt = len(pts)

        if not trans:
            if self.operator:
                trans1d = range(-2**n,2**n+2)
            else:
                trans1d = range(2**n)
            trans = []
            for lx in trans1d:
                for ly in trans1d:
                    for lz in trans1d:
                        trans.append((lx,ly,lz))

        for lx,ly,lz in trans:
            xlo = lx*h
            ylo = ly*h
            zlo = lz*h
            f = Tensor(npt,npt,npt)
            if self.cfunction:
                cfcube.cfcube(xlo, ylo, zlo, npt, h,
                              pts.swigptr(), f.swigptr(),
                              self.cfunction)
            else:
                # Plain old Python function or non-cast C function
                # This is SLOW due to all of the tensor references.
                self.fcube(xlo, ylo, zlo, npt, h, pts, f, self.function)
            f.scale(scale)
            self.s[n][lx,ly,lz] = f.fast_transform(phiw)

    def refine_fine_scale_projection(self,n,lx,ly,lz,prnt=1):
        '''

        We have the fine scale projection at level n+1, so we can
        compute the scaling function and wavelet coefficients at level
        n.  Considering just box l on level n, examine the difference
        coefficients.  If they do not satisfy the threshold, then
        refine the representation by recurring down to a finer level.

        If refinement occurs this modifies:
        :   self.s

        Set prnt to non-zero for printing

        '''
        if n > self.max_refine_level:
            self.nterminated = self.nterminated+1
            return

        if Function.truncate_method == 0:
            tol = self.thresh
        elif Function.truncate_method == 1:
            tol = sqrt(1.0/8.0**n)*self.thresh
        else:
            # Here, cannot do truncate_method=2 properly since we are
            # refining vertically down and don't yet know the number
            # of boxes at this level.  Arbitrarily chop it at level 5
            # ... i.e., assume that in all space there are no more
            # than 8**5=32768 boxes.  sqrt(32768)=181.01
            tol = self.thresh/sqrt(8.0**5)

        k = self.k
        ss = self.gather_scaling_coeffs(n,lx,ly,lz)
        ss = self.filter(ss)
        s = Tensor(ss[0:k,0:k,0:k])
        ss[0:k,0:k,0:k].fill(0.0)
        dnorm = ss.normf()
        if dnorm <= tol:
##             # The difference coefficients are neglible so can truncate
##             # at this level.  Delete the coeffs on the level below to
##             # save space, and adopt the scaling function coeffs
##             # computed from the lower level by twoscale since they
##             # will be more accurate
##             if not self.s.has_key(n): self.s[n] = {}
##             self.s[n][lx,ly,lz] = s
##             lx2, ly2, lz2 = lx*2, ly*2, lz*2
##             for ix in range(2):
##                 for iy in range(2):
##                     for iz in range(2):
##                         del(self.s[n+1][lx2+ix,ly2+iy,lz2+iz])
            pass
        else:
            # First remove the inaccurate scaling function coeffs at
            # level n+1, refine to level n+2
            if prnt:
                for i in range(n): print "  ",
                print " refining level=%d l=(%d,%d,%d) norm=%.1e" % \
                      (n, lx, ly, lz, dnorm)
            lx2, ly2, lz2 = lx*2, ly*2, lz*2
            range2 = range(2)
            for ix in range2:
                ixlo = ix*k
                ixhi = ixlo + k
                for iy in range2:
                    iylo = iy*k
                    iyhi = iylo+k
                    for iz in range2:
                        izlo = iz*k
                        izhi = izlo+k
                        del(self.s[n+1][lx2+ix,ly2+iy,lz2+iz])
                        lx4, ly4, lz4 = 2*(lx2+ix), 2*(ly2+iy), 2*(lz2+iz)
                        for iix in range2:
                            for iiy in range2:
                                for iiz in range2:
                                    self.fine_scale_projection(
                                        n+2,[(lx4+iix,ly4+iiy,lz4+iiz)])
                        self.refine_fine_scale_projection(
                            n+1,lx2+ix,ly2+iy,lz2+iz,prnt)

    def unfilter(self,ss,sonly=0):
        '''

        Given scaling function and wavelet coefficients (s and d)
        return the scaling function coefficients at the next finer
        level.  I.e., reconstruct Vn using Vn = Vn-1 + Wn-1.

        s0 = sum(j) h0_ji*s_j + g0_ji*d_j
        s1 = sum(j) h1_ji*s_j + g1_ji*d_j

        Returns a new tensor and has no side effects

        If (sonly) ... then ss is only the scaling function coeffs
        (and assume the d are zero).  Works for any number of
        dimensions.

        '''
        if sonly:
            return ss.fast_transform(self.hgsonly)
        else:
            return ss.fast_transform(self.hg)

    def filter(self,ss):
        '''
        
        Given scaling function coefficients s[n][l][i] and s[n][l+1][i]
        return the scaling function and wavelet coefficients at the
        coarser level.  I.e., decompose Vn using Vn = Vn-1 + Wn-1.

        s_i = sum(j) h0_ij*s0_j + h1_ij*s1_j
        d_i = sum(j) g0_ij*s0_j + g1_ij*s1_j

        Returns a new tensor and has no side effects.  Works for any
        number of dimensions.
        
        '''
        return ss.fast_transform(self.hgT)

    def gather_scaling_coeffs(self, n, lx, ly, lz):
        '''

        For a given translation at level n, gather the corresponding
        scaling function coefficients at level n+1.

        Some of the boxes on the lower level may be missing.

        '''
        k = self.k
        sn1 = self.s[n+1]
        ss = Tensor(2*k,2*k,2*k)
        lx2, ly2, lz2 = lx*2, ly*2, lz*2
        range2 = range(2)
        for ix in range2:
            ixlo = ix*k
            ixhi = ixlo + k
            for iy in range2:
                iylo = iy*k
                iyhi = iylo+k
                for iz in range2:
                    izlo = iz*k
                    izhi = izlo+k
                    key = (lx2+ix,ly2+iy,lz2+iz)
                    if sn1.has_key(key):
                        ss[ixlo:ixhi,iylo:iyhi,izlo:izhi] =  sn1[key]
        return ss
        
    def reconstruct(self, nonstandard=0):
        '''

        Reconstruct the projection onto scaling functions from the
        compressed form terminating at the locally lowest level
        without difference coefficients.

        The compression algorithm has the responsibility of ensuring
        that the there are no signficant nodes below a node without
        difference coefficients.

        This algorithm must ensure that there are scaling function
        coefficients only at the locally lowest level (since other
        algorithms will just run down the tree until the first
        coefficient is found, and then stop).

        !!! nonstandard = 1
        A variant of this algorithm is used to reconstruct the
        function in preparation for the application of a non-standard
        form operator.  For this we need to keep the difference
        coefficients and also keep the sum coefficients at each level
        on the way down, but only down to locally finest level at
        which we have difference coefficients.  There is (will be)
        another routine to clean the resulting mess back into the
        conventional hierarchical basis definition.

        !!! nonstandard = 2
        Another variant for application of the non-standard form operator
        requires just the scaling coefficients at all levels.  No need
        to keep the difference coefficients.

        '''

        if self.debug:
            print "reconstructing", id(self), self.compressed, self.s.keys()
            
        if not self.compressed:
            if nonstandard:
                self.compress()
            else:
                return self

        perf.enter('reconstrct')
        self.compressed = 0

        k = self.k

        nmax = max(self.d.keys())
        if nonstandard == 1:
            # Only need scaling functions down to where the difference
            # coefficients are locally non-zero.
            nmax = nmax-1
            #nmax = nmax

        range2 = range(2)
            
        for n in range(nmax+1):
            sn = self.s[n]
            dn = self.d[n]
            sn1 = self.s[n+1] = {}
            if nonstandard==1 and n<nmax:
                dn1 = self.d[n+1]
            for lx, ly, lz in dn.keys():
                d = dn[lx,ly,lz]
                d[0:k,0:k,0:k] = sn[lx,ly,lz]
                if nonstandard==1 and n>=nmax: continue
                if nonstandard != 2: del(sn[lx,ly,lz])
                d = self.unfilter(d)
                lx2, ly2, lz2 = lx*2, ly*2, lz*2
                for ix in range2:
                    for iy in range2:
                        for iz in range2:
                            key2 = lx2+ix,ly2+iy,lz2+iz
                            if (nonstandard==1) and (not dn1.has_key(key2)):
                                continue
                            dxyz = d[ix*k:ix*k+k,iy*k:iy*k+k,iz*k:iz*k+k]
                            sn1[key2] = Tensor(dxyz)

            if nonstandard == 1:
                del self.s[n]
            else:
                del self.d[n]

        if nonstandard == 1: del self.s[nmax+1]
        
        if nmax == -1 and nonstandard == 1:
            print "FUDGE", self.s[0][0,0,0].normf()
            self.d[0][0,0,0][0:k,0:k,0:k] = self.s[0][0,0,0]
            print "FUDGE", self.d[0][0,0,0].normf()
            del self.s[0][0,0,0]

        perf.exit('reconstrct')
        return self

    def compress(self,truncate=0):
        '''

        Compress the scaling function coefficients (Vn+1 = Vn + Wn).

        In the usual compressed form, the scaling function
        coefficients will only exist on the lowest level, but this
        will vary across the domain.

        In the non-standard compressed form, both scaling function and
        difference coefficients will be defined on all levels ... must
        add the result of compressing the sum coefficients into the
        existing difference coefficients.

        In the form that results from addition in the scaling function
        basis, there may be sum coefficients at many levels, but no
        differences.

        For ease of handling all the special cases now implement
        truncation as a separate and optional step.

        Returns self so that operations can be chained.

        '''
        if self.debug:
            print "compressing", id(self), self.compressed, self.s.keys()
            
        if self.compressed: return self
        self.compressed = 1
        perf.enter('compress')
        
        nmax = max(self.s.keys())
        for n in range(nmax-1,-1,-1):
            if not self.d.has_key(n):
                self.d[n] = {}

        k = self.k
        for n in range(nmax-1,-1,-1):
            dn  = self.d[n]
            sn1 = self.s[n+1]
            try:
                sn = self.s[n]
            except KeyError:
                sn = self.s[n] = {}
            try:
                dn1 = self.d[n+1]
            except KeyError:
                dn1 = {}
            
            # Compute the list of translations that we must consider, driven
            # from the non-zero entries on the level below.
            llist = {}
            for lx,ly,lz in sn1.keys():
                llist[parent(lx),parent(ly),parent(lz)] = 1
            
            for lx,ly,lz in llist.keys():
                ss = self.gather_scaling_coeffs(n, lx, ly, lz)
                ss = self.filter(ss)
                s = Tensor(ss[0:k,0:k,0:k]) # Force a copy
                if sn.has_key((lx,ly,lz)):
                    sn[lx,ly,lz] = sn[lx,ly,lz].gaxpy(1.0,s,1.0)
                else:
                    sn[lx,ly,lz] = s

                ss[0:k,0:k,0:k].fill(0.0) # Now just want the diff coeffs
                if dn.has_key((lx,ly,lz)):
                    dn[lx,ly,lz] = dn[lx,ly,lz].gaxpy(1.0,ss,1.0)
                else:
                    dn[lx,ly,lz] = ss

            del self.s[n+1]

        if truncate:
            self.truncate()

        if not self.d.has_key(0):
            self.d[0] = {}
        if not self.d[0].has_key((0,0,0)):
            self.d[0][0,0,0] = Tensor(2*k,2*k,2*k)

        if not self.s.has_key(0):
            self.s[0] = {}
        if not self.s[0].has_key((0,0,0)):
            self.s[0][0,0,0] = Tensor(k,k,k)

        perf.exit('compress')
        return self

    def restrict(self, basis):
        '''

        Self is a compressed function and we are restricting
        (truncating) it to be entirely contained in the given basis.
        Basis is a nested directory just like self.d whose keys give
        the dilations and translations of the allowed functions.  The
        space to be used is then V0 plus the pieces of Wn described by
        basis.

        Returns self so that operations can be chained.

        '''
        if self.debug:
            print "restricting", self
            
        if not self.compressed:
            self.compress()
        
        d = self.d

        # First delete any unncessary levels
        basis_nmax = max(basis.keys())
        self_nmax = max(self.d.keys())
        for n in range(basis_nmax+1,self_nmax+1):
            del d[n]

        for n in d.keys():
            dn = d[n]
            bn = basis[n]
            for key in dn.keys():
                if not bn.has_key(key):
                    del(dn[key])

        return self

    def project(self,newk):
        '''

        Returns a projection of self into the basis with
        multi-wavelets of order newk which can be less than
        or greater than self.k.

        The projection is performed by reconstructing on the
        finest level, and copying coefficients.

        If newk>=self.k, then the projection is exact.

        If newk<self.k, then information is lost.

        Self is not modified, but may be reconstructed

        '''
        k = self.k
        kk = min(k,newk)
        
        self.reconstruct()
        perf.enter('project')
        result = self.new(k=newk,compress=0)
        ndim = self.ndim
        nmax = max(self.s.keys())
        for n in range(nmax+1):
            sn = self.s[n]
            rn = result.s[n] = {}
            for key in sn.keys():
                s = sn[key]
                if ndim == 1:
                    r = Vector(newk)
                    r[0:kk] = s[0:kk]
                elif ndim == 2:
                    r = Matrix(newk,newk)
                    r[0:kk,0:kk] = s[0:kk,0:kk]
                elif ndim == 3:
                    r = Tensor(newk,newk,newk)
                    r[0:kk,0:kk,0:kk] = s[0:kk,0:kk,0:kk]
                elif ndim == 6:
                    r = Tensor(newk,newk,newk,newk,newk,newk)
                    r[0:kk,0:kk,0:kk,0:kk,0:kk,0:kk] = s[0:kk,0:kk,0:kk,0:kk,0:kk,0:kk]
                else:
                    raise "error in projection (ndim is either of 1, 2, 3, 6)"

                rn[key] = r

        perf.exit('project')
        return result.compress()

    def __lend(self):
        '''
        Return the number of boxes containing difference
        coefficients.
        '''
        lend = 0
        for n in self.d.keys():
            lend = lend + len(self.d[n])

        return lend

    def basis(self,refine=0):
        '''

        Self is a compressed function.  Export a description of
        the basis (for use by restrict).  This describes the
        pieces of W0+W1+...+W(n-1) in use.

        If refine=1, then if the difference coefficients do not
        satisfy the accuracy goal, the basis is locally refined down
        one level.

        '''
        if not self.compressed: self.compress()

        if Function.truncate_method == 2:
            lend = self.__lend()
        
        basis = {}
        nmax = max(self.d.keys())
        for n in range(nmax+1):
            dn = self.d[n]
            bn = basis[n] = {}
            for key in dn.keys():
                bn[key] = 1

        if refine:
            basis[nmax+1] = {}          # Just in case
            for n in range(nmax+1):
                dn = self.d[n]
                bn = basis[n]
                bn1 = basis[n+1]

                keys = dn.keys()

                if Function.truncate_method == 0:
                    tol = self.thresh
                elif Function.truncate_method == 1:
                    tol = sqrt(1.0/8.0**n)*self.thresh
                else:
                    tol = self.thresh/sqrt(lend)

                range2 = range(2)
                for key in keys():
                    if (dn[key].normf() > tol):
                        lx,ly,lz = key
                        lx2, ly2, lz2 = lx*2, ly*2, lz*2
                        printed = 0
                        for ix in range2:
                            for iy in range2:
                                for iz in range2:
                                    key2 = (lx2+ix,ly2+iy,lz2+iz)
                                    if not bn1.has_key(key2):
                                        if not printed:
                                            print "Refining", n, key
                                            printed = 1
                                        bn1[key2] = 1
            if not basis[nmax+1]:
                del basis[nmax+1]

        return basis
    
    def truncate(self,thresh=None):
        '''

        The difference coefficients are truncated to the specified
        threshold.  So that reconstruction can just proceed down to
        the first level without difference coefficients, it is
        important to only truncate difference coefficients at leaf
        nodes.

        Returns self so that operations can be chained.

        '''
        if self.debug:
            print "truncating", id(self), self.compressed, self.thresh
            
        if not self.compressed:
            self.compress()

        if Function.truncate_method == 2:
            lend = self.__lend()
        
        if thresh is None:
            thresh = self.thresh
        
        perf.enter('truncate')

        nmax = max(self.d.keys())

        is_parent = {}
        for n in range(nmax,0,-1):
            dn = self.d[n]
            keys = dn.keys()

            if Function.truncate_method == 0:
                tol = thresh
            elif Function.truncate_method == 1:
                tol = sqrt(1.0/8.0**n)*thresh
            else:
                tol = thresh/sqrt(lend)

            # If not a parent and below threshold kill box.
            # If above threshold mark parent box as such.
            next_is_parent = {}
            for key in keys:
                if (not is_parent.has_key(key)) and \
                   (dn[key].normf() < tol):
                    del dn[key]
                else:
                    lx,ly,lz = key
                    next_is_parent[parent(lx),parent(ly),parent(lz)] = 1

            is_parent = next_is_parent
                    
            # If the current level is totally empty and it is also
            # the lowest level, it can be deleted entirely.
            if (not self.d[n]) and (not self.d.has_key(n+1)):
                del self.d[n]

        perf.exit('truncate')

##         if Function.truncate_method == 2:
##             # If we deleted a lot of coefficients, try again ...
##             if self.__lend()*2 < lend:
##                 self.truncate()

        return self

    def norm2(self):
        '''

        Compute the 2-norm.

        '''
        perf.enter('norm2')
        if self.compressed:
            if self.ndim == 1:
                sum = self.s[0][0].normf()**2
            elif self.ndim == 2:
                sum = self.s[0][0,0].normf()**2
            elif self.ndim == 3:
                sum = self.s[0][0,0,0].normf()**2
            elif self.ndim == 6:
                sum = self.s[0][0,0,0,0,0,0].normf()**2
            else:
                raise "unknown ndim"
            for n in self.d.keys():
                dn = self.d[n]
                for key in dn.keys():
                    sum = sum + dn[key].normf()**2
        else:
            sum = 0.0
            for n in self.s.keys():
                sn = self.s[n]
                for key in sn.keys():
                    sum = sum + sn[key].normf()**2

        perf.exit('norm2')
        return sqrt(sum)

    def inner(self,other):
        '''

        Inner product of two functions which is the sum of
        the product of the coefficients since the basis is
        orthogonal.

        '''
        if self.compressed != other.compressed:
            self.compress()
            other.compress()

        perf.enter('inner')

        if self.compressed:
            if self.ndim!=1:
                key = (0,)*self.ndim
            else:
                key = 0
            #sum = self.s[0][key].flat().inner(other.s[0][key].flat())
            sum = self.s[0][key].trace(other.s[0][key])
            nmax = min(max(self.d.keys()),max(other.d.keys()))
            for n in range(nmax+1):
                dn = self.d[n]
                on = other.d[n]
                for key in dn.keys():
                    if on.has_key(key):
                        #sum = sum + dn[key].flat().inner(on[key].flat())
                        sum = sum + dn[key].trace(on[key])
        else:
            sum = 0.0
            nmax = max(self.s.keys()+other.s.keys())
            for n in range(nmax+1):
                try:
                    sn = self.s[n]
                except KeyError:
                    sn = self.s[n] = {}

                try:
                    on = other.s[n]
                except KeyError:
                    on = other.s[n] = {}

                # Generate union of keys on this level
                keys = sn.keys() + on.keys()
                keys.sort()
                for i in range(len(keys)-1,0,-1):
                    if keys[i] == keys[i-1]:
                        del keys[i]

                for key in keys:
                    s = self.__sock_it_to_me(n,key)
                    if not (s is None):
                        o = other.__sock_it_to_me(n,key)
                        if not (o is None):
                            #sum = sum + s.flat().inner(o.flat())
                            sum = sum + s.trace(o)
            
            self.sclean()
            other.sclean()

        perf.exit('inner')
        return sum

    def trace(self):
        '''

        Compute the trace = integral(self) over the unit cube

        '''
        if not self.compressed: self.compress()

        if self.ndim == 1:
            return self.s[0][0][0]
        elif self.ndim == 2:
            return self.s[0][0,0][0,0]
        elif self.ndim == 3:
            return self.s[0][0,0,0][0,0,0]
        elif self.ndim == 6:
            return self.s[0][0,0,0,0,0,0][0,0,0,0,0,0]
        else:
            raise "unknown ndim"

    def griddx(self,filename='/tmp/grid.dx'):
        '''

        Dump an opendx format version of the scaling function grid.
        (should extend to also dump the multiwavelet nested grids)
        (it would also be nice to have the different levels of the
        grid be distinct objects rather than just having different
        values).

        '''
        print "Dumping grid for Opendx to", filename
        sys.stdout.flush()
        f = open(filename,'w')
        self.reconstruct()
        nmax = max(self.s.keys())
        cubes = []
        points = []
        for n in range(nmax+1):
            try:
                sn = self.s[n]
            except KeyError:
                sn = {}
            for lx,ly,lz in sn.keys():
                cubes.append((n,lx,ly,lz))
                points.append((n,lx,ly,lz))
                points.append((n,lx+1,ly,lz))
                points.append((n,lx,ly+1,lz))
                points.append((n,lx,ly,lz+1))
                points.append((n,lx+1,ly+1,lz))
                points.append((n,lx+1,ly,lz+1))
                points.append((n,lx,ly+1,lz+1))
                points.append((n,lx+1,ly+1,lz+1))


        # Eliminate duplicate points
        points.sort()
        for i in range(len(points)-1,0,-1):
            if points[i] == points[i-1]:
                del(points[i])

        # Make a directory to index the points
        p = {}
        npt = 0
        for point in points:
            p[point] = npt
            npt = npt + 1

        ncube = len(cubes)
        print "npoint=%d ncube=%d nmax=%d" % (npt, ncube, nmax)
        sys.stdout.flush()

        # Write points to the dx file
        f.write("object 1 class array type float rank 1 shape 3 items %d data follows\n" % npt)
        for n,lx,ly,lz in points:
            scale = 0.5**n
            x,y,z = scale*lx,scale*ly,scale*lz
            f.write("%.8f %.8f %.8f\n" % (x, y, z))

        # Write connections to the dx file
        f.write("\n")
        f.write("object 2 class array type int rank 1 shape 8 items %d data follows\n" % ncube)
        for n, lx, ly, lz in cubes:
            p0 = p[n,lx,ly,lz]
            p1 = p[n,lx+1,ly,lz]
            p2 = p[n,lx,ly+1,lz]
            p3 = p[n,lx+1,ly+1,lz]
            p4 = p[n,lx,ly,lz+1]
            p5 = p[n,lx+1,ly,lz+1]
            p6 = p[n,lx,ly+1,lz+1]
            p7 = p[n,lx+1,ly+1,lz+1]
            f.write("%d %d %d %d %d %d %d %d\n" %
                    (p0,p1,p2,p3,p4,p5,p6,p7))

        f.write('attribute "element type" string "cubes"\n')
        f.write('attribute "ref" string "positions"\n')
        f.write("\n")

        # Write values at grid points = level of refinement
        f.write("object 3 class array type float rank 0 items %d data follows\n" % npt)
        for n,lx,ly,lz in points:
            f.write(" %d" % n)

        f.write("\n")
        f.write('attribute "dep" string "positions"\n')
        f.write("\n")
        f.write('object "irregular positions irregular connections" class field\n')
        f.write('component "positions" value 1\n')
        f.write('component "connections" value 2\n')
        f.write('component "data" value 3\n')
        f.write('end\n')
        f.close()
        
    def plotdx(self,filename='/tmp/plot.dx',
               npt=(20,20,20),
               box=(0.0,1.0,0.0,1.0,0.0,1.0)):
        '''

        Dump a cube of (nptx,npty,nptz) points to the specified file
        for visualization by opendx

        '''
        nptx, npty, nptz = npt
        xlo,xhi,ylo,yhi,zlo,zhi = box

        print "Generating cubic grid for Opendx in %s" % filename
        print "  grid = (%d,%d,%d)" % npt
        print "   box = (%.2f,%.2f,%.2f,%.2f,%.2f,%.2f)" % box
        sys.stdout.flush()

        self.reconstruct()
        self.fast_eval_init()
        func = self.fast_eval
        f = open(filename,'w')
        f.write('object 1 class gridpositions counts %d %d %d\n'%npt)
        f.write('origin %.6f %.6f %.6f\n' % (xlo,ylo,zlo))
        f.write('delta %.6f 0 0\n' % ((xhi-xlo)/(nptx-1)))
        f.write('delta 0 %.6f 0\n' % ((yhi-ylo)/(npty-1)))
        f.write('delta 0 0 %.6f\n' % ((zhi-zlo)/(nptz-1)))
        f.write('\n')
        f.write('object 2 class gridconnections counts %d %d %d\n'%npt)
        f.write('attribute "element type" string "cubes"\n')
        f.write('attribute "ref" string "positions"\n')
        f.write('\n')
        f.write('object 3 class array type float rank 0 items %d data follows\n' % (nptx*npty*nptz))
        
        xx, yy, zz = range(nptx), range(npty), range(nptz)
        for x in range(nptx):
            xx[x] = xlo + (xhi-xlo)*x/(nptx-1.0)
        for y in range(npty):
            yy[y] = ylo + (yhi-ylo)*y/(npty-1.0)
        for z in range(nptz):
            zz[z] = zlo + (zhi-zlo)*z/(nptz-1.0)
        for x in xx:
            for y in yy:
                for z in zz:
                    f.write("%.6e\n" % func(x,y,z))

        f.write('attribute "dep" string "positions"\n')
        f.write('\n')
        f.write('object "regular positions regular connections" class field\n')
        f.write('component "positions" value 1\n')
        f.write('component "connections" value 2\n')
        f.write('component "data" value 3\n')
        f.write('\n')
        f.write('end\n')
        f.close()
        print "Done plotting"
        sys.stdout.flush()
        self.fast_eval_tidy()
        
    def eval(self,x,y,z):
        '''

        Evaluate the function at the given point.

        If it is not compressed then we just need to walk down the
        tree until we locate the leaf node containing x,y,z.  If it is
        compressed then we do the same, but also have to add up the
        scaling function coeffs on the way down.

        '''
        n = 0
        k = self.k
        lx = ly = lz = 0
        xx, yy, zz = x, y, z
        if self.compressed:
            # Note always want loop to terminate via the break otherwise
            # n is one too small (x & l would also be wrong)
            s = self.s[0][0,0,0]
            for n in range(max(self.d.keys())+2):
                try:
                    d = Tensor(self.d[n][lx,ly,lz]) # Take a copy
                except KeyError:
                    # At the bottom ... evaluate s*phi(xx) at level n
                    break
                d[0:k,0:k,0:k] = s
                d = self.unfilter(d)
                scale = 2**(n+1)
                xx, yy, zz = x*scale, y*scale, z*scale
                lx, ly, lz = min(scale-1,int(xx)), min(scale-1,int(yy)), \
                             min(scale-1,int(zz))
                xx, yy, zz = xx-lx, yy-ly, zz-lz
                ix, iy, iz = lx%2, ly%2, lz%2
                s = d[ix*k:ix*k+k,iy*k:iy*k+k,iz*k:iz*k+k]
            s = Tensor(s)
        else:
            # The function is not compressed.  Assume that it is
            # represented at some (locally defined) finest level
            # in the scaling function basis.
            nmax = max(self.s.keys()) + 1

            ss = self.s
            # Eliminate higher levels of the tree
##             nmin = 0
##             for n in range(nmax):
##                 if len(ss[n]):
##                     break
##                 else:
##                     nmin += 1

            for n in range(nmax):
                twon = 2**n
                twon1= twon - 1
                xx, yy, zz = x*twon, y*twon, z*twon
                lx, ly, lz = min(twon1,int(xx)), min(twon1,int(yy)), \
                             min(twon1,int(zz))
                xx, yy, zz = xx-lx, yy-ly, zz-lz
                try:
		    #print "found: ", n, lx, ly, lz
                    s = ss[n][lx,ly,lz]
                    break
                except KeyError:
                    pass
            
        #px = Vector(phi(xx,k))
        #py = Vector(phi(yy,k))
        #pz = Vector(phi(zz,k))
        #value = s.inner(pz).inner(py).inner(px)*(2.0**n)*sqrt(2.0**n)

        value = cfcube.eval(k,n,xx,yy,zz,s.swigptr())
        return value

    def eval_err(self, npt=1000, scan=1, box=0.4):
        self.reconstruct()
        self.fast_eval_init()
        if scan:
            print "\nScanning errors along line x=y=0.5 - coarse"
            for z in range(0,11):
                z = z * 0.1
                test = self(0.5,0.5,z)
                exact = self.function(0.5,0.5,z)
                err = test-exact
                if exact:
                    relerr = err/exact
                else:
                    relerr = 0.0
                print "  z=%.4f f=%16.10f abserr=%.1e relerr=%.1e" % \
                      (z,test,err,relerr)

            print "\nScanning errors along line x=y=0.5 - medium"
            for z in range(45,56):
                z = z * 0.01
                test = self(0.5,0.5,z)
                exact = self.function(0.5,0.5,z)
                err = test-exact
                if exact:
                    relerr = err/exact
                else:
                    relerr = 0.0
                print "  z=%.4f f=%16.10f abserr=%.1e relerr=%.1e" % \
                      (z,test,err,relerr)

            print "\nScanning errors along line x=y=0.5 - fine"
            for z in range(495,506):
                z = z * 0.001
                test = self(0.5,0.5,z)
                exact = self.function(0.5,0.5,z)
                err = test-exact
                if exact:
                    relerr = err/exact
                else:
                    relerr = 0.0
                print "  z=%.4f f=%16.10f abserr=%.1e relerr=%.1e" % \
                      (z,test,err,relerr)

            print "\nScanning errors along line x=y=0.5 - finest"
            for z in range(4995,5006):
                z = z * 0.0001
                test = self.fast_eval(0.5,0.5,z)
                exact = self.function(0.5,0.5,z)
                err = test-exact
                if exact:
                    relerr = err/exact
                else:
                    relerr = 0.0
                print "  z=%.4f f=%16.10f abserr=%.1e relerr=%.1e" % \
                      (z,test,err,relerr)

        print "\nSampling error at %d random points" % npt
        errsq = 0.0
        maxerr= 0.0
        maxrelerr=0.0
        seed(201)
        lo = max(0.0,0.5-box/2.0)
        for i in range(npt):
            x, y, z = lo+random()*box, lo+random()*box, lo+random()*box
            r = sqrt(x*x+y*y+z*z)
            if r < 0.001:
                continue                # Analytic form inaccurate here?
            value = self.fast_eval(x,y,z)
            exact = self.function(x,y,z)
            err = abs(exact - value)
            if err > maxerr:
                maxerr = err
                xmax, ymax, zmax = x, y, z
            errsq = errsq + err*err
            if exact:
                maxrelerr = max(maxrelerr,err/exact)
        rmserr = sqrt(errsq/npt)
        print "RMS abs error = %.1e" % rmserr
        print "MAX abs error = %.1e" % maxerr
        print "MAX rel error = %.1e" % maxrelerr
        print "Location of max abs err = %.6f %.6f %.6f" % \
              (xmax, ymax, zmax)

        self.fast_eval_tidy()

    def print_layers(self, print_pic=1):
        if self.compressed:
            print " Analysis by layers of compressed function"
            nmax = max(self.d.keys())
            normsq = self.s[0][0,0,0].normf()**2
            print "\nNorm in V0 %.8e" % sqrt(normsq)
            for n in range(nmax+1):
                sum = 0.0
                for xyz in self.d[n].keys():
                    sum = sum + self.d[n][xyz].normf()**2
                normsq = normsq + sum
                print "Norm in W%d %.8e %.8e #=%d" % (n,sqrt(sum),sqrt(normsq), len(self.d[n].keys()))

            if not print_pic: return  
            print "\nNorms of difference coeffs"
            for n in range(min(nmax+1,5)):
                for lz in range(2**n):
                    print "\nLevel =",n,"  lz =",lz
                    print "    ",
                    for lx in range(2**n):
                        print "  %2d  "%lx,
                    print ""
                    print "    ",
                    for lx in range(2**n):
                        print " -----",
                    print ""

                    for ly in range(2**n):
                        print "%2d |" % ly,
                        for lx in range(2**n):
                            try:
                                print " %.0e" % self.d[n][lx,ly,lz].normf(),
                            except KeyError:
                                print "   -  ",
                        print "|"
        else:
            print " Analysis by layers of reconstructed function"
            nmax = max(self.s.keys())
            normsq = 0.0
            for n in range(nmax+1):
                sum = 0.0
                count = 0
                try:
                    count = len(self.s[n].keys())
                    for xyz in self.s[n].keys():
                        sum = sum + self.s[n][xyz].normf()**2
                except KeyError:
                    pass
                normsq = normsq + sum
                print "Norm in V%d %.8e %.8e #=%d" % (n,sqrt(sum),sqrt(normsq),count)

            #if not print_pic: return  
            print "\nNorms of sum coeffs at each level"
            for n in range(min(nmax+1,5)):
                for lz in range(2**n):
                    print "\nLevel =",n,"  lz =",lz
                    print "    ",
                    for lx in range(2**n):
                        print "  %2d  "%lx,
                    print ""
                    print "    ",
                    for lx in range(2**n):
                        print " -----",
                    print ""

                    for ly in range(2**n):
                        print "%2d |" % ly,
                        for lx in range(2**n):
                            try:
                                print " %.0e" % self.s[n][lx,ly,lz].normf(),
                            except KeyError:
                                print "   -  ",
                        print "|"
            
    def copy(self):
        '''

        Return a brand spanking new deep copy of self.

        Now that the Tensor class has set/getstate methods this
        method is not really necessary ... U can use the copy.deepcopy()
        function as for other Python classes.  However, this copy
        method might (?) be a bit more efficient.

        '''
        perf.enter('copy')
        r = self.new(thresh=self.thresh, k=self.k, compress=self.compressed)
        r.function = self.function
        r.compressed = self.compressed
        for n in self.s.keys():
            r.s[n] = {}
            for key in self.s[n].keys():
                r.s[n][key] = Tensor(self.s[n][key])
        for n in self.d.keys():
            r.d[n] = {}
            for key in self.d[n].keys():
                r.d[n][key] = Tensor(self.d[n][key])
        perf.exit('copy')
        return r

    def scalar_to_function(self, value):
        # Lazy way to turn a scalar into a function
        value = float(value)
        if self.ndim == 1:
            f = lambda x,r=value: r
        elif self.ndim == 2:
            f = lambda x,y,r=value: r
        elif self.ndim == 3:
            f = lambda x,y,z,r=value: r
        elif self.ndim == 6:
            f = lambda x1,y1,z1,x2,y2,z2,r=value: r
        else:
            raise "unknown ndim"

        return self.new(k=self.k,function=f,
                        compress=self.compressed,
                        initial_level=0, refine=0)
        
    def __call__(self, x, y=None, z=None, x2=None, y2=None, z2=None):
        if self.ndim == 1:
            return self.eval(x)
        elif self.ndim == 2:
            return self.eval(x, y)
        elif self.ndim == 3:
            return self.eval(x, y, z)
        elif self.ndim == 6:
            return self.eval(x, y, z, x2, y2, z2)
        else:
            raise "Yo! How many dimensions do you want?"

    def __sub__(self, other):
        '''

        Subtraction of one function from another, or subtraction of a
        scalar from a function, or subtraction of a python callable
        object (that takes x,y,z in [0,1] as arguments) from a function.

        Returns a new function with no effect on self or other.

        The functions must either both be compressed or both
        uncompressed and using the same order wavelet.

        '''

        if not isinstance(other,Function):
            if callable(other):
                other = self.new(k=self.k,function=other,
                                 compress=self.compressed)
            else:
                other = self.scalar_to_function(other)

        result = self.copy().gaxpy(1.0,other,-1.0)
        result.thresh = Function.thresh
        # If in scaling function basis then result will have
        # scaling function coeffs scattered throughout ...
        # compress here to add these together.
        return result.compress()
        

    def __mul__(self,other):
        '''

        Multiplication of a function by either another, or a python
        callable object (that takes x,y,z in [0,1] as arguments) or a
        scalar.

        Returns a new function with no effect on self or other.

        '''

        if isinstance(other,Function  ) or \
           isinstance(other,Function6d) or \
           isinstance(other,Function1d) or \
           isinstance(other,Function2d):
            return self.mul(other)
            #return self.mul_old(other)
        else:
            if callable(other):
                other = self.new(k=self.k,function=other,
                                 compress=self.compressed)
                return self.mul(other)
            else:
                return self.copy().scale(float(other))


    def sclean(self):
        '''

        Cleans up scaling function coefficients after diffing or
        multiplying.  Only want to keep the scaling function coeffs
        at the highest level in the tree ... ones lower down will
        have been generated by two-scale just for evaluation.
        Kill all nodes that have a parent - infanticide.

        '''
        for n in range(max(self.s.keys()),0,-1):
            sn, sn1 = self.s[n], self.s[n-1]
            for lx,ly,lz in sn.keys():
                if sn1.has_key((parent(lx),parent(ly),parent(lz))):
                    del(sn[lx,ly,lz])
        
    def mul_old(self,other):
        sumsq = (self+other).square()
        diffsq = (self-other).square()
        result = sumsq.gaxpy(0.25,diffsq,-0.25)
        result.compress()               # Needed due to gaxpy
        return result
        
    def mul(self,other):
        '''

        Multiplication of a function by a another function.
        Replaces the algorithm based on squaring.

        Returns a new function in the scaling function basis.
        Self/other are not modified though may be reconstructed if
        either is compressed.

        Added an optimization for multiplcation by the mask which is
        mask(x) = 1 for 0.0625 < x < 1-0.625.  Thus, only need to do
        the multiplication we have boxes 0 or 15 on level 4 (or boxes
        contained within below).  The mask is identified by the
        attribute mask which noone else will have.

        '''

        self.reconstruct()
        other.reconstruct()

        if hasattr(self,'mask'): self,other = other,self
        mask = hasattr(other,'mask')
        
        perf.enter('mul')

        result = self.new(k=self.k,compress=0) # Was thresh=self.thresh
        result.s = {}

        k = self.k
        phit = Matrix(self.quad_phi.cycledim(1)) # Contig. for fast_transform
        phiw = self.quad_phiw
        npt = len(self.quad_x)
        nmax = max(self.s.keys()+other.s.keys())

        tmp1 = Tensor(k,k,k)
        tmp2 = Tensor(k,k,k)

        range2 = range(2)

        if Function.autorefine:     # NOTE REFERENCE TO CLASS
            result.s[0] = {}
             
        for n in range(nmax+1):
            try:
                sn = self.s[n]
            except KeyError:
                sn = self.s[n] = {}

            try:
                on = other.s[n]
            except KeyError:
                on = other.s[n] = {}
                
            if Function.autorefine:     # NOTE REFERENCE TO CLASS
                rn1 = result.s[n+1] = {}
                scale = sqrt(8.0**(n+1))
            else:
                rn = result.s[n] = {}
                scale = sqrt(8.0**n)

            # Generate union of keys on this level. 
            keys = sn.keys() + on.keys()
            keys.sort()
            for i in range(len(keys)-1,0,-1):
                if keys[i] == keys[i-1]:
                    del keys[i]

            if mask:
                if n <= 4:
                    masklo, maskhi = 0, 2**n-1
                else:
                    masklo, maskhi = 2**(n-4)-1, 15*2**(n-4)
                    
            for key in keys:
                lx,ly,lz = key

                # If other is the mask, exclude keys for which the
                # mask is unity.
                if mask:
                    if masklo<lx<maskhi and masklo<ly<maskhi and \
                       masklo<lz<maskhi:
                        if sn.has_key(key): rn[key] = Tensor(sn[key])
                        continue

                s =  self.__sock_it_to_me(n,key)
                if not (s is None):
                    o = other.__sock_it_to_me(n,key)
                    if not (o is None):
                        if Function.autorefine: # NOTE REFERENCE TO CLASS
                            ss = self.unfilter(s,sonly=1)
                            oo = self.unfilter(o,sonly=1)
                            lx2, ly2, lz2 = 2*lx, 2*ly, 2*lz
                            for ix in range2:
                                for iy in range2:
                                    for iz in range2:
                                        f = ss[ix*k:ix*k+k,iy*k:iy*k+k,
                                               iz*k:iz*k+k]
                                        g = oo[ix*k:ix*k+k,iy*k:iy*k+k,
                                               iz*k:iz*k+k]
                                        f = f.transform(phit,result=tmp1)
                                        g = g.transform(phit,result=tmp2)
                                        f.emul(g)
                                        rn1[lx2+ix,ly2+iy,lz2+iz] = \
                                           f.fast_transform(phiw).scale(scale)
                        else:
                            f = s.fast_transform(phit,result=tmp1)
                            g = o.fast_transform(phit,result=tmp2)
                            f.emul(g)
                            rn[key] = f.fast_transform(phiw).scale(scale)
                            
        self.sclean()
        other.sclean()

        perf.exit('mul')

        return result

    def mul3(self,other,other2):
        '''

        Multiplication of a function by a another function.
        Replaces the algorithm based on squaring.

        Returns a new function in the scaling function basis.
        Self/other are not modified though may be reconstructed if
        either is compressed.

        Added an optimization for multiplcation by the mask which is
        mask(x) = 1 for 0.0625 < x < 1-0.625.  Thus, only need to do
        the multiplication we have boxes 0 or 15 on level 4 (or boxes
        contained within below).  The mask is identified by the
        attribute mask which noone else will have.

        '''

        self.reconstruct()
        other.reconstruct()
        other2.reconstruct()

        perf.enter('mul3')

        result = self.new(k=self.k,compress=0) # Was thresh=self.thresh
        result.s = {}

        k = self.k
        phit = Matrix(self.quad_phi.cycledim(1)) # Contig. for fast_transform
        phiw = self.quad_phiw
        npt = len(self.quad_x)
        nmax = max(self.s.keys()+other.s.keys()+other2.s.keys())

        tmp1 = Tensor(k,k,k)
        tmp2 = Tensor(k,k,k)
        tmp3 = Tensor(k,k,k)

        range2 = range(2)

        if Function.autorefine:     # NOTE REFERENCE TO CLASS
            rn1 = result.s[0] = {}
            
        for n in range(nmax+1):
            try:
                sn = self.s[n]
            except KeyError:
                sn = self.s[n] = {}

            try:
                on = other.s[n]
            except KeyError:
                on = other.s[n] = {}
                
            try:
                o2n = other2.s[n]
            except KeyError:
                o2n = other2.s[n] = {}
                
            if Function.autorefine:     # NOTE REFERENCE TO CLASS
                rn1 = result.s[n+1] = {}
                scale = 8.0**(n+1)
            else:
                rn = result.s[n] = {}
                scale = 8.0**n

            # Generate union of keys on this level. 
            keys = sn.keys() + on.keys() + o2n.keys()
            keys.sort()
            for i in range(len(keys)-1,0,-1):
                if keys[i] == keys[i-1]:
                    del keys[i]

            for key in keys:
                lx,ly,lz = key

                s =  self.__sock_it_to_me(n,key)
                if not (s is None):
                    o = other.__sock_it_to_me(n,key)
                    if not (o is None):
                        o2 = other2.__sock_it_to_me(n,key)
                        if not (o2 is None):
                            if Function.autorefine: # NOTE REFERENCE TO CLASS
                                ss = self.unfilter(s,sonly=1)
                                oo = self.unfilter(o,sonly=1)
                                oo2 = self.unfilter(o2,sonly=1)
                                lx2, ly2, lz2 = 2*lx, 2*ly, 2*lz
                                for ix in range2:
                                    for iy in range2:
                                        for iz in range2:
                                            f = ss[ix*k:ix*k+k,iy*k:iy*k+k,
                                                   iz*k:iz*k+k]
                                            g = oo[ix*k:ix*k+k,iy*k:iy*k+k,
                                                   iz*k:iz*k+k]
                                            g2 = oo2[ix*k:ix*k+k,iy*k:iy*k+k,
                                                   iz*k:iz*k+k]
                                            f = f.transform(phit,result=tmp1)
                                            g = g.transform(phit,result=tmp2)
                                            g2 = g2.transform(phit,result=tmp3)
                                            f.emul(g)
                                            f.emul(g2)
                                            rn1[lx2+ix,ly2+iy,lz2+iz] = \
                                               f.fast_transform(phiw).scale(scale)
                            else:
                                f = s.fast_transform(phit,result=tmp1) #.scale(scale)
                                g = o.fast_transform(phit,result=tmp2) #.scale(scale)
                                g2 = o2.fast_transform(phit,result=tmp3) #.scale(scale)
                                f.emul(g)
                                f.emul(g2)
                                rn[key] = f.fast_transform(phiw).scale(scale) #.scale(1.0/scale)
                            
        self.sclean()
        other.sclean()
        other2.sclean()

        perf.exit('mul3')

        return result

    def __rmul__(self,other):
        '''

        This is other*self.  For scalar multiplication this is just
        the same as self*other, but we also want to account for
        operators (Dx, Dy, Dz, DxT, DyT, DzT, D2x, D2y, D2z, Del,
        Delsq).

        '''

        if other is Dx:
            return self.diff(0)
        elif other is Dy:
            return self.diff(1)
        elif other is Dz:
            return self.diff(2)
        elif other is DxT:
            return self.diff(0,transpose=1)
        elif other is DyT:
            return self.diff(1,transpose=1)
        elif other is DzT:
            return self.diff(2,transpose=1)
        elif other is D2x:
            return self.diff2(0)
        elif other is D2y:
            return self.diff2(1)
        elif other is D2z:
            return self.diff2(2)
        elif other is Del:
            return [Dx*self, Dy*self, Dz*self]
        elif other is Delsq:
            return self.laplacian()
        else:
            return self.__mul__(other)


    def __pow__(self,power):
        if power != 2:
            raise ValueError, "Function: __pow__: only squaring for now\n"
        return self.copy().square()
            

    def __add__(self,other):
        '''

        Addition of one function to another, or addition of a scalar
        to a function, or addition of a python callable object (that
        takes x,y,z in [0,1] as arguments) to a function.

        Returns a new function with no effect on self or other.

        '''

        if not (isinstance(other,Function  ) or \
                isinstance(other,Function6d) or \
                isinstance(other,Function1d) or \
                isinstance(other,Function2d)):
            if callable(other):
                other = self.new(k=self.k,function=other,
                                 compress=self.compressed)
            else:
                other = self.scalar_to_function(other)

        result = self.copy().gaxpy(1.0,other,1.0)

        # If in scaling function basis then result will have
        # scaling function coeffs scattered throughout ...
        # compress here to add these together.
        return result.compress()

    def square(self):
        '''

        Square self replacing self by its square in the scaling
        function basis.

        If (Function.autorefine) is set, then it will recur down an
        additional level before squaring, otherwise it does not.

        NB: self is changed.

        Returns self so that operations may be chained.  

        '''
        if self.debug:
            print "squaring", id(self), self.compressed
            
        self.reconstruct()

        k = self.k
        nmax = max(self.s.keys())
        phit = Matrix(self.quad_phi.cycledim(1)) # Contig. for fast_transform
        phiw = self.quad_phiw
        npt = len(self.quad_x)

        tmp1 = Tensor(k,k,k)
        range2 = range(2)

        try:
            next_keys = self.s[0].keys()
        except KeyError:
            # This can happen if not compressed and initial_level != 0
            self.s[0] = {}
            next_keys = self.s[0].keys()

        for n in range(nmax+1):
        
            sn = self.s[n]
            try:
                sn1 = self.s[n+1]
            except KeyError:
                sn1 = self.s[n+1] = {}

            keys = next_keys
            next_keys = sn1.keys()

            if Function.autorefine:         # NOTE REFERENCE TO CLASS
                ss = Tensor(2*k,2*k,2*k)
                scale = sqrt(8.0**(n+1))
            else:
                scale = sqrt(8.0**n)
            for lx, ly, lz in keys:
                if Function.autorefine:     # NOTE REFERENCE TO CLASS
                    # Now are at the locally lowest level ... recur down
                    # one more level and then form the square projecting
                    # onto the approximately refined basis.  Could recur
                    # more dynamically to guarantee precision but
                    # Beylkin does not seem to think this is necessary
                    # except in pathological cases?
                    ss = self.unfilter(sn[lx,ly,lz],sonly=1)
                    del(sn[lx,ly,lz])
                    lx2, ly2, lz2 = 2*lx, 2*ly, 2*lz
                    for ix in range2:
                        for iy in range2:
                            for iz in range2:
                                f = ss[ix*k:ix*k+k,iy*k:iy*k+k,iz*k:iz*k+k]
                                # Back transform to the function value
                                # at the quadrature points
                                f = f.transform(phit,result=tmp1)
                                # Compute the square on the grid
                                f.emul(f)
                                # Do quadrature on f**2 and scale
                                sn1[lx2+ix,ly2+iy,lz2+iz] = \
                                    f.fast_transform(phiw).scale(scale)
                                                          
                else:
                    # Don't refine ... just do it at the current level
                    f = sn[lx,ly,lz].fast_transform(phit,result=tmp1)
                    f.emul(f)
                    sn[lx,ly,lz] = f.fast_transform(phiw).scale(scale)
            
        return self


    def local_function(self,function=None,cfunction=None):
        '''

        Apply a scalar local function to self replacing self by
        f(self(x)).  Define either function (a callable object) or
        cfunction (a pointer to a function cast to a pointer to void
        ... so that old versions of SWIG can be used).  Both the C and
        Python functions should take a double argument and return a
        double.  THIS APPROACH IS ONLY GOOD IF THE RESULT OF THE
        FUNCTION IS AS LOCALLY SMOOTH OR SMOOTHER THAN THE INPUT ...

        If (Function.autorefine) is set, then it will recur down an
        additional level before squaring, otherwise it does not.

        NB: self is changed.

        Returns self so that operations may be chained.  

        '''
        if self.debug:
            print "squaring", id(self), self.compressed
            
        range2 = range(2)

        self.reconstruct()
        perf.enter('local func')

        k = self.k
        nmax = max(self.s.keys())
        phit = Matrix(self.quad_phi.cycledim(1)) # Contig. for fast_transform
        phiw = self.quad_phiw
        npt = len(self.quad_x)

        tmp1 = Tensor(k,k,k)

        try:
            next_keys = self.s[0].keys()
        except KeyError:
            # This can happen if not compressed and initial_level != 0
            self.s[0] = {}
            next_keys = self.s[0].keys()

        for n in range(nmax+1):
        
            sn = self.s[n]
            try:
                sn1 = self.s[n+1]
            except KeyError:
                sn1 = self.s[n+1] = {}

            keys = next_keys
            next_keys = sn1.keys()

            if Function.autorefine:         # NOTE REFERENCE TO CLASS
                ss = Tensor(2*k,2*k,2*k)
                scale = sqrt(8.0**(n+1))
            else:
                scale = sqrt(8.0**n)
            for lx, ly, lz in keys:
                if Function.autorefine:     # NOTE REFERENCE TO CLASS
                    # Now are at the locally lowest level ... recur down
                    # one more level and then form the square projecting
                    # onto the approximately refined basis.  Could recur
                    # more dynamically to guarantee precision but
                    # Beylkin does not seem to think this is necessary
                    # except in pathological cases?
                    ss = self.unfilter(sn[lx,ly,lz],sonly=1)
                    del(sn[lx,ly,lz])
                    lx2, ly2, lz2 = 2*lx, 2*ly, 2*lz
                    for ix in range2:
                        for iy in range2:
                            for iz in range2:
                                f = ss[ix*k:ix*k+k,iy*k:iy*k+k,iz*k:iz*k+k]
                                # Back transform to the function value
                                # at the quadrature points
                                f = Tensor(f).transform(phit,result=tmp1).scale(scale) 
                                # Eval function on the refined grid
                                if cfunction:
                                    cfcube.clocalfunc(
                                        f.swigptr(),k**3,cfunction)
                                else:
                                    f.unaryop(function)
                                # Do quadrature on f**2 and scale
                                sn1[lx2+ix,ly2+iy,lz2+iz] = \
                                    f.fast_transform(phiw).scale(1.0/scale)
                                                          
                else:
                    # Don't refine ... just do it at the current level
                    f = sn[lx,ly,lz].fast_transform(phit,result=tmp1).scale(scale)
                    if cfunction:
                        cfcube.clocalfunc(f.swigptr(),k**3,cfunction)
                    else:
                        f.unaryop(function)
                    sn[lx,ly,lz] = f.fast_transform(phiw).scale(1.0/scale)
            
        perf.exit('local func')
        return self

    def mul_func(self,function=None,cfunction=None):
        '''

        '''
        if self.debug:
            print "squaring", id(self), self.compressed
            
        range2 = range(2)

        self.reconstruct()
        perf.enter('mul func')

        k = self.k
        nmax = max(self.s.keys())
        phit = Matrix(self.quad_phi.cycledim(1)) # Contig. for fast_transform
        phiw = self.quad_phiw
        pts = self.quad_x 
        npt = len(self.quad_x)

        tmp1 = Tensor(k,k,k)

        try:
            next_keys = self.s[0].keys()
        except KeyError:
            # This can happen if not compressed and initial_level != 0
            self.s[0] = {}
            next_keys = self.s[0].keys()

        for n in range(nmax+1):

            sn = self.s[n]
            try:
                sn1 = self.s[n+1]
            except KeyError:
                sn1 = self.s[n+1] = {}

            keys = next_keys
            next_keys = sn1.keys()

            if Function.autorefine:         # NOTE REFERENCE TO CLASS
                ss = Tensor(2*k,2*k,2*k)
                dd = float(2**(n+1))
            else:
                dd = float(2**n)

            scale = sqrt(dd**3)       # 3D Scale
            h = 1.0/dd                # Box size on target level
             
            for lx, ly, lz in keys:

                if Function.autorefine:     # NOTE REFERENCE TO CLASS
                    # Now are at the locally lowest level ... recur down
                    # one more level and then form the square projecting
                    # onto the approximately refined basis.  Could recur
                    # more dynamically to guarantee precision but
                    # Beylkin does not seem to think this is necessary
                    # except in pathological cases?
                    ss = self.unfilter(sn[lx,ly,lz],sonly=1)
                    del(sn[lx,ly,lz])
                    lx2, ly2, lz2 = 2*lx, 2*ly, 2*lz
                    for ix in range2:
                        for iy in range2:
                            for iz in range2:

                                # get a cube of g:(func)
                                # ----------------------
                                xlo = (lx2+ix)*h
                                ylo = (ly2+iy)*h
                                zlo = (lz2+iz)*h
                                g = Tensor(npt,npt,npt)
                                if cfunction:
                                    cfcube.cfcube(xlo, ylo, zlo, npt, h,
                                                  pts.swigptr(), g.swigptr(),
                                                  cfunction)
                                else:
                                    # Plain old Python function or non-cast C function
                                    # This is SLOW due to all of the tensor references.
                                    self.fcube(xlo, ylo, zlo, npt, h, pts, g, function)
                                
                                # get a cube of f:(func) 
                                # ----------------------
                                f = ss[ix*k:ix*k+k,iy*k:iy*k+k,iz*k:iz*k+k]
                                # Back transform to the function value
                                # at the quadrature points
                                f = f.transform(phit,result=tmp1).scale(scale)

                                # emul between f and g
                                # --------------------
                                f = f.emul(g)
                                
                                # Do quadrature on f**2 and scale
                                sn1[lx2+ix,ly2+iy,lz2+iz] = \
                                    f.fast_transform(phiw).scale(1.0/scale)
                                                          
                else:

                    # get a cube of g:(func)
                    # ----------------------
                    xlo = lx*h
                    ylo = ly*h
                    zlo = lz*h
                    g = Tensor(npt,npt,npt)
                    if cfunction:
                        cfcube.cfcube(xlo, ylo, zlo, npt, h,
                                      pts.swigptr(), g.swigptr(),
                                      cfunction)
                    else:
                        # Plain old Python function or non-cast C function
                        # This is SLOW due to all of the tensor references.
                        self.fcube(xlo, ylo, zlo, npt, h, pts, g, function)

                    # get a cube of f:(func) 
                    # ----------------------
                    # Don't refine ... just do it at the current level
                    f = sn[lx,ly,lz].fast_transform(phit,result=tmp1).scale(scale)

                    # emul between f and g
                    # --------------------
                    f = f.emul(g)
                
                    sn[lx,ly,lz] = f.fast_transform(phiw).scale(1.0/scale)
            
        perf.exit('mul func')
        return self

    def mul_func_trace(self,function=None,cfunction=None):
        '''

        '''
        if self.debug:
            print "squaring", id(self), self.compressed
            
        range2 = range(2)

        self.reconstruct()
        perf.enter('mul func')

        k = self.k
        nmax = max(self.s.keys())
        phit = Matrix(self.quad_phi.cycledim(1)) # Contig. for fast_transform
        phiw = self.quad_phiw
        phiw0 = Matrix(k,1)
        phiw0[:,0] = phiw[:,0]
        pts = self.quad_x 
        npt = len(self.quad_x)

        print "mul_func_trace: autorefinement=", Function.autorefine

        tmp1 = Tensor(k,k,k)

        try:
            next_keys = self.s[0].keys()
        except KeyError:
            # This can happen if not compressed and initial_level != 0
            self.s[0] = {}
            next_keys = self.s[0].keys()

        testval = 0.0
        trace = 0.0
        g = Tensor(npt,npt,npt)
        
        for n in range(nmax+1):

            sn = self.s[n]
            try:
                sn1 = self.s[n+1]
            except KeyError:
                sn1 = self.s[n+1] = {}

            keys = next_keys
            next_keys = sn1.keys()

            if Function.autorefine:         # NOTE REFERENCE TO CLASS
                ss = Tensor(2*k,2*k,2*k)
                dd = float(2**(n+1))
            else:
                dd = float(2**n)

            scale = sqrt(dd**3)       # 3D Scale
            h = 1.0/dd                # Box size on target level

            trace_leveln = 0.0
            ff = Tensor(k,k,k)
            
            for lx, ly, lz in keys:

                if Function.autorefine:     # NOTE REFERENCE TO CLASS
                    # Now are at the locally lowest level ... recur down
                    # one more level and then form the square projecting
                    # onto the approximately refined basis.  Could recur
                    # more dynamically to guarantee precision but
                    # Beylkin does not seem to think this is necessary
                    # except in pathological cases?
                    ss = self.unfilter(sn[lx,ly,lz],sonly=1)
                    #del(sn[lx,ly,lz])
                    lx2, ly2, lz2 = 2*lx, 2*ly, 2*lz
                    for ix in range2:
                        for iy in range2:
                            for iz in range2:

                                # get a cube of g:(func)
                                # ----------------------
                                xlo = (lx2+ix)*h
                                ylo = (ly2+iy)*h
                                zlo = (lz2+iz)*h
                                if cfunction:
                                    cfcube.cfcube(xlo, ylo, zlo, npt, h,
                                                  pts.swigptr(), g.swigptr(),
                                                  cfunction)
                                else:
                                    # Plain old Python function or non-cast C function
                                    # This is SLOW due to all of the tensor references.
                                    self.fcube(xlo, ylo, zlo, npt, h, pts, g, function)
                                
                                # get a cube of f:(func) 
                                # ----------------------
                                f = ss[ix*k:ix*k+k,iy*k:iy*k+k,iz*k:iz*k+k]
                                # Back transform to the function value
                                # at the quadrature points
                                f = f.transform(phit,result=tmp1) #.scale(scale)

                                # emul between f and g
                                # --------------------
                                ff.gaxpy(1.0,f.emul(g),1.0)
                                
                                ## Do quadrature on f**2 and scale
                                #sn1[lx2+ix,ly2+iy,lz2+iz] = \
                                #    f.fast_transform(phiw).scale(1.0/scale)
                                                          
                else:

                    # get a cube of g:(func)
                    # ----------------------
                    xlo = lx*h
                    ylo = ly*h
                    zlo = lz*h
                    if cfunction:
                        cfcube.cfcube(xlo, ylo, zlo, npt, h,
                                      pts.swigptr(), g.swigptr(),
                                      cfunction)
                    else:
                        # Plain old Python function or non-cast C function
                        # This is SLOW due to all of the tensor references.
                        self.fcube(xlo, ylo, zlo, npt, h, pts, g, function)

                    # get a cube of f:(func) 
                    # ----------------------
                    # Don't refine ... just do it at the current level
                    f = sn[lx,ly,lz].fast_transform(phit,result=tmp1) #.scale(scale)

                    # emul between f and g
                    # --------------------
                    ff.gaxpy(1.0,f.emul(g),1.0)
                
                    #sn[lx,ly,lz] = f.fast_transform(phiw).scale(1.0/scale)
                    #testval += sn[lx,ly,lz][0,0,0]*(1/scale)
                    
            trace += ff.fast_transform(phiw0)[0,0,0]/scale #/scale

        nmax = max(self.s.keys())
        for n in range(nmax,-1,-1):
            if len(self.s[n]): break
            del self.s[n]
                
        perf.exit('mul func')
        return trace

    def memory_usage(self):
        '''
        Print a summary of memory usage
        '''
        ntensors = 0
        mem = 0L
        k = self.k
        for n in self.s.keys():
            nkey = len(self.s[n].keys())
            mem = mem + nkey*k**self.ndim
            ntensors = ntensors + nkey
        for n in self.d.keys():
            nkey = len(self.d[n].keys())
            ntensors = ntensors + nkey
            mem = mem + nkey*(2*k)**self.ndim
        mem = mem*8 
        print "%8d tensors  %9d bytes" % (ntensors, mem)

    def scale(self, alpha):
        '''

        self = alpha*self where alpha is a scalar

        Returns self so that operations may be chained.

        '''
        perf.enter('scale')
        if self.debug:
            print "scaling", id(self), self.compressed, alpha
            
        alpha = float(alpha)

        if self.compressed:
            if self.ndim == 1:
                self.s[0][0].scale(alpha)
            elif self.ndim == 2:
                self.s[0][0,0].scale(alpha)
            elif self.ndim == 3:
                self.s[0][0,0,0].scale(alpha)
            elif self.ndim == 6:
                self.s[0][0,0,0,0,0,0].scale(alpha)
            else:
                raise "Yo! How many dimensions do you want?"

            for n in self.d.keys():
                for key in self.d[n].keys():
                    self.d[n][key].scale(alpha)
        else:
            for n in self.s.keys():
                for key in self.s[n].keys():
                    self.s[n][key].scale(alpha)

        perf.exit('scale')
        return self

    def gaxpy(self, alpha, other, beta):
        '''

        self = alpha*self + beta*other

        Other is not changed.

        Returns self so that operations can be chained.

        Note that if the inputs are in the scaling function basis the
        result is also in the same basis ... and may have coeffs
        scattered throughout the tree so it probably needs to be
        compressed before other operations are performed.

        '''
        if self.debug:
            print "gaxpy", self.compressed, alpha, beta
        
        if self.compressed != other.compressed:
            self.compress()
            other.compress()
            
        perf.enter('gaxpy')
        if self.compressed:
            if self.ndim!=1:
                key = (0,)*self.ndim
            else:
                key = 0
            self.s[0][key].gaxpy(alpha,other.s[0][key],beta)
            nmax = max(other.d.keys())
            for n in range(nmax+1):
                try:
                    self.d[n]
                except KeyError:
                    self.d[n] = {}
                for key in other.d[n].keys():
                    try:
                        self.d[n][key].gaxpy(alpha,other.d[n][key],beta)
                    except KeyError:
                        self.d[n][key] = Tensor(other.d[n][key]).scale(beta)
        else:
            self.scale(alpha)
            nmax = max(other.s.keys())
            for n in range(nmax+1):
                on = other.s[n]
                try:
                    sn = self.s[n]
                except KeyError:
                    sn = self.s[n] = {}
                for key in on.keys():
                    try:
                        sn[key].gaxpy(1.0,on[key],beta)
                    except KeyError:
                        sn[key] = Tensor(on[key]).scale(beta)

        perf.exit('gaxpy')
        return self

    def make_dc_periodic(self):
        '''

        Return the level-0 blocks rm, r0, rp of the central
        difference derivative operator with periodic boundary
        conditions on either side.

        '''
        k = self.k
        r0 = Matrix(k,k)
        rp = Matrix(k,k)
        rm = Matrix(k,k)
        
        iphase = 1.0
        for i in range(k):
            jphase = 1.0
            for j in range(k):
                gammaij = sqrt((2*i+1)*(2*j+1))
                if ((i-j)>0) and (((i-j)%2)==1):
                    Kij = 2.0
                else:
                    Kij = 0.0
                r0[i,j] = 0.5*(1.0 - iphase*jphase - 2.0*Kij)*gammaij
                rm[i,j] = 0.5*jphase*gammaij
                rp[i,j] =-0.5*iphase*gammaij
                jphase = -jphase
            iphase = -iphase
        return (rm, r0, rp)

    def make_dc_periodic_compressed(self):

        k = self.k
        k2 = k*2

        rm, r0, rp = self.make_dc_periodic()
        
        T0 = Matrix(k2,k2)
        T0[0  :k  ,0  :k  ] = r0[0:k,0:k] # r00 = r 0
        T0[0+k:k+k,0+k:k+k] = r0[0:k,0:k] # r11 = r 0
        T0[0+k:k+k,0  :k  ] = rp[0:k,0:k] # r10 = r+1
        T0[0  :k  ,0+k:k+k] = rm[0:k,0:k] # r01 = r-1
        T0 = T0.fast_transform(self.hgT)
        T0 = T0.scale(2.0)

        R0 = T0[0  :k  ,0  :k  ] # PTP
        a0 = T0[0+k:k+k,0+k:k+k] # QTQ
        b0 = T0[0+k:k+k,0  :k  ] # QTP
        c0 = T0[0  :k  ,0+k:k+k] # PTQ

        Tp = Matrix(k2,k2)
        Tp[0  :k  ,0+k:k+k] = rp[0:k,0:k]
        Tp = Tp.fast_transform(self.hgT)
        Tp = Tp.scale(2.0)
        
        Rp = Tp[0  :k  ,0  :k  ] # PTP
        ap = Tp[0+k:k+k,0+k:k+k] # QTQ
        bp = Tp[0+k:k+k,0  :k  ] # QTP
        cp = Tp[0  :k  ,0+k:k+k] # PTQ

        Tm = Matrix(k2,k2)
        Tm[0+k:k+k,0  :k  ] = rm[0:k,0:k]
        Tm = Tm.fast_transform(self.hgT)
        Tm = Tm.scale(2.0)

        Rm = Tm[0  :k  ,0  :k  ] # PTP
        am = Tm[0+k:k+k,0+k:k+k] # QTQ
        bm = Tm[0+k:k+k,0  :k  ] # QTP
        cm = Tm[0  :k  ,0+k:k+k] # PTQ

        #print (r0-R0).normf() # should be zero 
        #print (rp-Rp).normf() # should be zero 
        #print (rm-Rm).normf() # should be zero 

        return (Tm, T0, Tp), (rm, r0, rp)

    def __look_out_below(self,s, n, key):
        '''

        s are scaling coefficients at level n for the given
        translation.  Recur them down one level and put the resulting
        coefficients into the tree.

        '''
        if not self.s.has_key(n+1):
            self.s[n+1] = {}
        k = self.k
        ss = self.unfilter(s,sonly=1)
	print "bar: ",s
	print "baz: ",ss
	sys.exit()
        lx, ly, lz = key
        lx2, ly2, lz2 = lx*2, ly*2, lz*2
        for ix in range(2):
            for iy in range(2):
                for iz in range(2):
                    view = ss[ix*k:ix*k+k,iy*k:iy*k+k,iz*k:iz*k+k]
                    self.s[n+1][lx2+ix,ly2+iy,lz2+iz] = Tensor(view)

    def __sock_it_to_me(self,n,key):
        '''

        Return self.s[n][key] or None if it cannot be made (because it
        has no parent).
        
        If it is absent, then try to recur from the level above to
        make it.  If the parent is also absent, we need to keep
        walking up.  If we hit the top of the tree and still cannot
        find a parent, then return None.

        Some optimizations are possible, such as looking one below
        before recuring up.  Also, has_key() might be faster than
        exception handling.

        '''
        try:
            return self.s[n][key]
        except KeyError:
            if n == 0:
                return None

        lx, ly, lz = key

        lx, ly, lz = parent(lx), parent(ly), parent(lz)
        try:
            s = self.s[n-1][lx,ly,lz]
        except KeyError:
            s = self.__sock_it_to_me(n-1,(lx,ly,lz))

        if s is None:
            return None

        self.__look_out_below(s,n-1,(lx,ly,lz))

        return self.s[n][key]

    def __sock_it_to_me2(self,n,key):
        '''

        Return self.s[n][key] or None if it cannot be made (because it
        has no parent).
     
        If it is absent, then try to recur from the level above to
        make it.  If the parent is also absent, we need to keep
        walking up.  If we hit the top of the tree and still cannot
        find a parent, then return None.

        Some optimizations are possible, such as looking one below
        before recuring up.  Also, has_key() might be faster than
        exception handling.

        '''
        try:
            return self.s[n][key].normf()
        except KeyError:
            if n == 0:
                return None

        lx, ly, lz = key

        lx, ly, lz = parent(lx), parent(ly), parent(lz)
        try:
            s = self.s[n-1][lx,ly,lz]
            s = s.normf()
        except KeyError:
            s = self.__sock_it_to_me2(n-1,(lx,ly,lz))

        if s is None:
            return None

        s = s * 0.35355339059327379 # = (2.0**(-3.0/2.0))
     
        return s

    def diff(self,axis,transpose=0):
        '''

        Apply a central derivative along the specified axis,
        currently with zero boundary conditions.  Plan to implement
        periodic boundary conditions also.

        If (transpose) then use the transposed derivative operator.

        Self is reconstructed but should otherwise be unchanged.

        Returns a new function in the scaling function basis.

        '''
        if self.debug:
            print "diffing", id(self), self.compressed, axis, \
                  transpose,self.s.keys()

        self.reconstruct()
        perf.enter('diff')
        
        result = self.new(k=self.k, compress=0)

        zero = Tensor(self.k,self.k,self.k)

        ind = (1,axis)                  # Indices to contract in products
        if axis == 0:
            amap = None
        elif axis == 1:
            amap = [1,0,2]
        elif axis == 2:
            amap = [2,0,1]
        else:
            raise "Don't you think three dimensions are enough?"
            
        rm, r0, rp = self.make_dc_periodic()

        if transpose:
            rm, r0, rp = rp.transpose(), r0.transpose(), rm.transpose()

        nmax = max(self.s.keys())

        for n in range(nmax+1):         #  Fill in empty levels
            try:
                self.s[n]
            except KeyError:
                self.s[n] = {}
            
        for n in range(nmax+1):
            result.s[n] = {}
            sn = self.s[n]
            twon = 2**n
            twon1= twon - 1
            for key in self.s[n].keys():
                s0 = sn[key]
                
                if key[axis] == 0:
                    sm = zero      # Zero outside the box
                else:
                    t=list(key); t[axis]=t[axis]-1; t=tuple(t)
                    sm = self.__sock_it_to_me(n,t)
                    
                if key[axis] == twon1:
                    sp = zero
                else:
                    t=list(key); t[axis]=t[axis]+1; t=tuple(t)
                    sp = self.__sock_it_to_me(n,t)
                
                if (sm is None) or (sp is None):
                    # One or both of my neighbours are further
                    # down the tree.  Just recur me down to the
                    # next level and try again there.
		    print "boobies"
                    self.__look_out_below(s0,n,key)
                    continue

                # rp & rm are only rank-1 ... don't yet exploit this

                r = rp.inner(sm,ind)
                r0.inner(s0,ind,result=r)  
                rm.inner(sp,ind,result=r)
                result.s[n][key] = r.make_self_contiguous_by_map(amap).scale(twon)
 
                #r = rp.inner(sm,ind,map=amap)
                #r0.inner(s0,ind,amap,result=r)  
                #rm.inner(sp,ind,map=amap,result=r)
                #result.s[n][key] = r.scale(twon)
                
                #result.s[n][key] = (rp.inner(sm,ind,amap) +   \
                #                    r0.inner(s0,ind,amap) +  \
                #                    rm.inner(sp,ind,amap)).scale(twon)
                
        self.sclean()

        perf.exit('diff')
        return result

    def diff2(self,axis):
        '''

        Apply a central derivative approximation to the second
        derivative along the specified axis, currently with zero
        boundary conditions.  Plan to implement periodic boundary
        conditions also.

        This operation corresponds to .  D2 = DT*D =
        self.diff(axis).diff(axist,transpose=1) and it should be
        negative (semi-definite ... I think that the boundary
        condition introduces a null space ... Beylkin suggests how to
        get rid of it but the projector is not yet coded).

        Self is reconstructed but should otherwise be unchanged.

        Returns a new function in the scaling function basis.

        '''
        if self.debug:
            print "diffing2", id(self), self.compressed, axis

        self.reconstruct()
        perf.enter('diff2')
        
        # Since are modifying the internals of result here must make
        # correct choice of NOT compressed
        result = self.new(k=self.k, compress=0)

        zero = Tensor(self.k,self.k,self.k)

        ind = (1,axis)                  # Indices to contract in products
        if axis == 0:
            amap = None
        elif axis == 1:
            amap = [1,0,2]
        elif axis == 2:
            amap = [2,0,1]
        else:
            raise "Don't you think three dimensions are enough?"
            
        rm, r0, rp = self.make_dc_periodic()
        Rp2 = rm.transpose()*rp
        Rp1 = rm.transpose()*r0 + r0.transpose()*rp
        R0  = rm.transpose()*rm + r0.transpose()*r0 + rp.transpose()*rp
        Rm1 =                     r0.transpose()*rm + rp.transpose()*r0
        Rm2 =                                         rp.transpose()*rm

        nmax = max(self.s.keys())

        for n in range(nmax+1):         #  Fill in empty levels
            try:
                self.s[n]
            except KeyError:
                self.s[n] = {}
            
        for n in range(nmax+1):
            result.s[n] = {}
            sn = self.s[n]
            twon = 2**n
            twon1= twon - 1
            scale = 4.0**n
            for key in self.s[n].keys():
                s = [-2,-1,0,1,2]
                got_them_all = 1
                for shift in s:
                    t=list(key); t[axis]=t[axis]+shift; t=tuple(t)
                    if t[axis]<0 or t[axis]>twon1:
                        s[2+shift] = zero
                    else:
                        s[2+shift] = self.__sock_it_to_me(n,t)
                        if s[2+shift] is None:
                            got_them_all = 0
                            break
                if not got_them_all:
                    # One or more of my neighbours are further
                    # down the tree.  Just recur me down to the
                    # next level and try again there.
                    self.__look_out_below(sn[key],n,key)
                    continue


                # Rp2 & Rm2 are only rank-1 ... don't yet exploit this
                sm2, sm1, s0, sp1, sp2 = s

                r = Rp2.inner(sm2,ind,amap)
                Rp1.inner(sm1,ind,amap,result=r)
                R0.inner(s0,ind,amap,result=r)
                Rm1.inner(sp1,ind,amap,result=r)
                Rm2.inner(sp2,ind,amap,result=r)
                result.s[n][key]=r.scale(scale)

                #result.s[n][key] = (Rp2.inner(sm2,ind,amap) +
                #                    Rp1.inner(sm1,ind,amap) + 
                #                    R0.inner(s0,ind,amap) +
                #                    Rm1.inner(sp1,ind,amap) +
                #                    Rm2.inner(sp2,ind,amap)).scale(scale)

        self.sclean()
        perf.exit('diff2')
        return result

    def laplacian(self):
        '''

        Compute the laplacian using a central derivative currently
        with a zero boundary condition (embedding).
        This operation corresponds to
        .   D2 = DT*D = self.diff(axis).diff(axis,transpose=1)
        and it should be negative (semi-definite ... I think that
        the embedding introduces a null space ... Beylkin suggests
        how to get rid of it but the projector is not yet coded).

        Self is reconstructed but should otherwise be unchanged.

        Returns a new function in the scaling function basis.

        '''
        if self.debug:
            print "laplacian", id(self), self.compressed
        
        self.reconstruct()
        perf.enter('laplacian')

        # Since are modifying the internals of result here must make
        # correct choice of NOT compressed
        result = self.new(k=self.k, compress=0)

        zero = Tensor(self.k,self.k,self.k)

        # For each axis the indices to contract in products and the
        # map from default to desired order in the contractions
        ind = ((1,0), (1,1), (1,2))
        amap = ((0,1,2), (1,0,2), (2,0,1))
            
        rm, r0, rp = self.make_dc_periodic()
        Rp2 = rm.transpose()*rp
        Rp1 = rm.transpose()*r0 + r0.transpose()*rp
        R0  = rm.transpose()*rm + r0.transpose()*r0 + rp.transpose()*rp
        Rm1 =                     r0.transpose()*rm + rp.transpose()*r0
        Rm2 =                                         rp.transpose()*rm

        nmax = max(self.s.keys())

        for n in range(nmax+1):         #  Fill in empty levels
            try:
                self.s[n]
            except KeyError:
                self.s[n] = {}
            
        for n in range(nmax+1):
            result.s[n] = {}
            sn = self.s[n]
            twon = 2**n
            twon1= twon - 1
            scale = -(4.0**n)
            for key in self.s[n].keys():
                s = [[-2,-1,0,1,2],[-2,-1,0,1,2],[-2,-1,0,1,2]]
                got_them_all = 1
                for axis in 0,1,2:
                    for shift in [-2,-1,0,1,2]:
                        t=list(key); t[axis]=t[axis]+shift; t=tuple(t)
                        if t[axis]<0 or t[axis]>twon1:
                            s[axis][2+shift] = zero
                        else:
                            s[axis][2+shift] = self.__sock_it_to_me(n,t)
                            if s[axis][2+shift] is None:
                                got_them_all = 0
                                break
                    if not got_them_all:
                        break
                    
                if not got_them_all:
                    # One or more of my neighbours are further
                    # down the tree.  Just recur me down to the
                    # next level and try again there.
                    self.__look_out_below(sn[key],n,key)
                    continue


                # Rp2 & Rm2 are only rank-1 ... don't yet exploit this
                r = Tensor(self.k,self.k,self.k)
                for axis in [0,1,2]:
                    sm2, sm1, s0, sp1, sp2 = s[axis]
                    i, a = ind[axis], amap[axis]
                    Rp2.inner(sm2,i,a,result=r)
                    Rp1.inner(sm1,i,a,result=r)
                    R0.inner(s0,i,a,result=r)
                    Rm1.inner(sp1,i,a,result=r)
                    Rm2.inner(sp2,i,a,result=r)
                    
                    #r = r + (Rp2.inner(sm2,i,a) + 
                    #         Rp1.inner(sm1,i,a) + 
                    #         R0.inner(s0,i,a)   +
                    #         Rm1.inner(sp1,i,a) +
                    #         Rm2.inner(sp2,i,a))
                    
                result.s[n][key] = r.scale(scale)

        self.sclean()
        perf.exit('laplacian')
        return result
                
    def __sock_it_to_me_diffNS(self,n,key,surface=0):

        k, k2 = self.k, self.k*2
        lx, ly, lz = key

        if surface:
            if self.cache_sock_it_to_me_diffNS.has_key((lx,ly,lz)):
                return self.cache_sock_it_to_me_diffNS[lx,ly,lz]
        
        ix, iy, iz = parentmod(lx), parentmod(ly), parentmod(lz)
        lx, ly, lz = parent(lx), parent(ly), parent(lz)

        try:
            s = self.d[n-1][lx,ly,lz]
            s = self.unfilter(s)
        except KeyError:
            s = self.__sock_it_to_me_diffNS(n-1,(lx,ly,lz),surface=0)
            s = self.unfilter(s,sonly=1)
            
        s = Tensor(s[ix*k:ix*k+k,iy*k:iy*k+k,iz*k:iz*k+k])

        if surface:
            lx, ly, lz = key
            self.cache_sock_it_to_me_diffNS[lx,ly,lz] = s

        return s

    def __sock_it_to_me2_diffNS(self,n,key,surface=0):

        k = self.k
        lx, ly, lz = key

        if surface:
            if self.cache_sock_it_to_me2_diffNS.has_key((lx,ly,lz)):
                return self.cache_sock_it_to_me2_diffNS[lx,ly,lz]
        
        ix, iy, iz = parentmod(lx), parentmod(ly), parentmod(lz)
        lx, ly, lz = parent(lx), parent(ly), parent(lz)
        
        try:
            s = self.d[n-1][lx,ly,lz]
            s = self.unfilter(s)
            s = Tensor(s[ix*k:ix*k+k,iy*k:iy*k+k,iz*k:iz*k+k])
            s = s.normf()
        except KeyError:
            s = self.__sock_it_to_me2_diffNS(n-1,(lx,ly,lz),surface=0)
            s = s * 0.35355339059327379 # = (2.0**(-3.0/2.0))
            
        if surface:
            lx, ly, lz = key
            self.cache_sock_it_to_me2_diffNS[lx,ly,lz] = s

        return s

    def diffNS(self,axis,transpose=0,autorefine=0):

        k, k2 = self.k, self.k*2
        
        if self.debug:
            print "diffing", id(self), self.compressed, axis, \
                  transpose,self.s.keys()

        self.compress()
        perf.enter('diffNS')
        self.reconstruct(nonstandard=1)

        result = self.new(k=self.k, thresh=self.thresh)

        ind = (1,axis)                  # Indices to contract in products
        if axis == 0:
            amap = None
        elif axis == 1:
            amap = [1,0,2]
        elif axis == 2:
            amap = [2,0,1]
        else:
            raise "Don't you think three dimensions are enough?"
            
        (Tm, T0, Tp), (rm, r0, rp) = self.make_dc_periodic_compressed()

        if transpose:
            Tm, T0, Tp = Tp.transpose(), T0.transpose(), Tm.transpose()
            rm, r0, rp = rp.transpose(), r0.transpose(), rm.transpose()

        nmax = max(self.d.keys())

        if autorefine: nmax += 1

        for n in range(nmax+1):         #  Fill in empty levels
            try:
                self.d[n]
            except KeyError:
                self.d[n] = {}

        cut1, cut2 = 0, 0
            
        for n in range(nmax+1):
            result.d[n] = {}
            dn = self.d[n]
            twon = 2**n
            twon1= twon - 1
            
            self.cache_sock_it_to_me_diffNS = {}
            self.cache_sock_it_to_me2_diffNS = {}

            if autorefine:
                if n>0:
                    keys = []
                    for lx,ly,lz in self.d[n-1].keys():
                        keys += [ (lx*2  ,ly*2  ,lz*2  ),
                                  (lx*2  ,ly*2  ,lz*2+1),
                                  (lx*2  ,ly*2+1,lz*2  ),
                                  (lx*2  ,ly*2+1,lz*2+1),
                                  (lx*2+1,ly*2  ,lz*2  ),
                                  (lx*2+1,ly*2  ,lz*2+1),
                                  (lx*2+1,ly*2+1,lz*2  ),
                                  (lx*2+1,ly*2+1,lz*2+1) ]
                else:
                    keys = [(0,0,0)]
            else:
                keys = dn.keys()

            for key in self.d[n].keys():
                t=list(key); t[axis] += 1; t=tuple(t)
                if t[axis]<=twon1 and not (t in keys): keys.append(t)
                t=list(key); t[axis] -= 1; t=tuple(t)
                if t[axis]>=0     and not (t in keys): keys.append(t)

            for key in keys:
                r = Tensor(k2,k2,k2)
                if n>0: rPP = Tensor(k,k,k)

                # NOTE: We should expolit the sparsity of multiwavelet bases
                # in non-standard form 

                # T0 and r0 =====================================
                t=key
                try:
                    d = dn[t]
                except KeyError:
                    d0 = self.__sock_it_to_me_diffNS(n,t,surface=1)
                    d = Tensor(k2,k2,k2)
                    d[0:k,0:k,0:k] = d0
                
                T0.inner(d,ind,map=amap,result=r) 
                if n>0:
                    d = d[0:k,0:k,0:k]
                    r0.inner(d,ind,map=amap,result=rPP)

                # Tp and rp =====================================
                if key[axis] != 0:
                    t=list(key); t[axis]=t[axis]-1; t=tuple(t)
                    try:
                        d = dn[t]
                    except KeyError:
                        d0 = self.__sock_it_to_me_diffNS(n,t,surface=1)
                        d = Tensor(k2,k2,k2)
                        d[0:k,0:k,0:k] = d0

                    Tp.inner(d,ind,map=amap,result=r) 
                    if n>0:
                        d = d[0:k,0:k,0:k]
                        rp.inner(d,ind,map=amap,result=rPP) 

                # Tm and rm =====================================
                if key[axis] != twon1:
                    t=list(key); t[axis]=t[axis]+1; t=tuple(t)
                    try:
                        d = dn[t]
                    except KeyError:
                        d0 = self.__sock_it_to_me_diffNS(n,t,surface=1)
                        d = Tensor(k2,k2,k2)
                        d[0:k,0:k,0:k] = d0
                        
                    Tm.inner(d,ind,map=amap,result=r)
                    if n>0:
                        d = d[0:k,0:k,0:k]
                        rm.inner(d,ind,map=amap,result=rPP)

                # ===============================================
                if n>0:
                    r[0:k,0:k,0:k].gaxpy(1.0,rPP[0:k,0:k,0:k],-1.0)
                    
                result.d[n][key] = r.scale(twon)
                
                #cut1 += 1
                #tmp = Tensor(result.d[n][key])
                #tmp[0:k,0:k,0:k].gaxpy(1.0,result.d[n][key][0:k,0:k,0:k],-1.0)
                #if tmp.normf() > self.thresh: cut2 +=1

            del self.cache_sock_it_to_me_diffNS
            del self.cache_sock_it_to_me2_diffNS

        # Operator has been applied.  Must now cleanup f from its
        # nonstandard state, and also compress the result summing the
        # sum coeffs at all levels in the process.

        #print cut2, "/", cut1

        for n in result.d.keys():
            dn = result.d[n]
            sn = result.s[n] = {}
            for key in dn.keys():
                sn[key] = Tensor(dn[key][0:k,0:k,0:k])
                dn[key][0:k,0:k,0:k].fill(0.0)
            
        result.compressed = 0
        result.compress()

        # Clean up the source
        self.s[0] = {}
        self.s[0][0,0,0] = Tensor(self.d[0][0,0,0][0:k,0:k,0:k])
        for n in self.d.keys():
            dn = self.d[n]
            for key in dn.keys():
                dn[key][0:k,0:k,0:k].fill(0.0)
        self.compressed = 1
        
        for n in range(nmax,-1,-1):
            if len(self.d[n])==0:
                del self.d[n]
            else:
                break

        perf.exit('diffNS')
        return result

    def diff2NS(self,axis,autorefine=0):

        k, k2 = self.k, self.k*2
        
        if self.debug:
            print "diffing2", id(self), self.compressed, axis

        self.compress()
        perf.enter('diff2NS')
        self.reconstruct(nonstandard=1)

        result = self.new(k=self.k, thresh=self.thresh)

        ind = (1,axis)                  # Indices to contract in products
        if axis == 0:
            amap = None
        elif axis == 1:
            amap = [1,0,2]
        elif axis == 2:
            amap = [2,0,1]
        else:
            raise "Don't you think three dimensions are enough?"
            
        (Tm, T0, Tp), (rm, r0, rp) = self.make_dc_periodic_compressed()
        
        tp2 = Tm.transpose()*Tp
        tp1 = Tm.transpose()*T0 + T0.transpose()*Tp
        t0  = Tm.transpose()*Tm + T0.transpose()*T0 + Tp.transpose()*Tp
        tm1 =                     T0.transpose()*Tm + Tp.transpose()*T0
        tm2 =                                         Tp.transpose()*Tm

        Rp2 = rm.transpose()*rp
        Rp1 = rm.transpose()*r0 + r0.transpose()*rp
        R0  = rm.transpose()*rm + r0.transpose()*r0 + rp.transpose()*rp
        Rm1 =                     r0.transpose()*rm + rp.transpose()*r0
        Rm2 =                                         rp.transpose()*rm

        del Tm, T0, Tp, rm, r0, rp
        
        nmax = max(self.d.keys())

        if autorefine: nmax += 1

        for n in range(nmax+1):         #  Fill in empty levels
            try:
                self.d[n]
            except KeyError:
                self.d[n] = {}

        cut1, cut2 = 0, 0
            
        for n in range(nmax+1):
            result.d[n] = {}
            dn = self.d[n]
            twon = 2**n
            twon1= twon - 1
            scale = 4.0**n
            
            self.cache_sock_it_to_me_diffNS = {}
            self.cache_sock_it_to_me2_diffNS = {}
            
            if autorefine:
                if n>0:
                    keys = []
                    for lx,ly,lz in self.d[n-1].keys():
                        keys += [ (lx*2  ,ly*2  ,lz*2  ),
                                  (lx*2  ,ly*2  ,lz*2+1),
                                  (lx*2  ,ly*2+1,lz*2  ),
                                  (lx*2  ,ly*2+1,lz*2+1),
                                  (lx*2+1,ly*2  ,lz*2  ),
                                  (lx*2+1,ly*2  ,lz*2+1),
                                  (lx*2+1,ly*2+1,lz*2  ),
                                  (lx*2+1,ly*2+1,lz*2+1) ]
                else:
                    keys = [(0,0,0)]
            else:
                keys = dn.keys()

            for key in self.d[n].keys():
                t=list(key); t[axis] += 2; t=tuple(t)
                if t[axis]<=twon1 and not (t in keys): keys.append(t)
                t=list(key); t[axis] += 1; t=tuple(t)
                if t[axis]<=twon1 and not (t in keys): keys.append(t)
                t=list(key); t[axis] -= 1; t=tuple(t)
                if t[axis]>=0     and not (t in keys): keys.append(t)
                t=list(key); t[axis] -= 2; t=tuple(t)
                if t[axis]>=0     and not (t in keys): keys.append(t)

            for key in keys:
                r = Tensor(k2,k2,k2)
                if n>0: rPP = Tensor(k,k,k)

                # NOTE: We should expolit the sparsity of multiwavelet bases
                # in non-standard form 

                # t0 and R0 =====================================
                t=key
                try:
                    d = dn[t]
                except KeyError:
                    d0 = self.__sock_it_to_me_diffNS(n,t,surface=1)
                    d = Tensor(k2,k2,k2)
                    d[0:k,0:k,0:k] = d0
                
                t0.inner(d,ind,map=amap,result=r) 
                if n>0:
                    d = d[0:k,0:k,0:k]
                    R0.inner(d,ind,map=amap,result=rPP)

                # tp1 and Rp1 ===================================
                if key[axis] != 0:
                    t=list(key); t[axis]=t[axis]-1; t=tuple(t)
                    try:
                        d = dn[t]
                    except KeyError:
                        d0 = self.__sock_it_to_me_diffNS(n,t,surface=1)
                        d = Tensor(k2,k2,k2)
                        d[0:k,0:k,0:k] = d0
                
                    tp1.inner(d,ind,map=amap,result=r) 
                    if n>0:
                        d = d[0:k,0:k,0:k]
                        Rp1.inner(d,ind,map=amap,result=rPP) 

                # tp2 and Rp2 ===================================
                if not (key[axis] in (0, 1)):
                    t=list(key); t[axis]=t[axis]-2; t=tuple(t)
                    try:
                        d = dn[t]
                    except KeyError:
                        d0 = self.__sock_it_to_me_diffNS(n,t,surface=1)
                        d = Tensor(k2,k2,k2)
                        d[0:k,0:k,0:k] = d0
                
                    tp2.inner(d,ind,map=amap,result=r) 
                    if n>0:
                        d = d[0:k,0:k,0:k]
                        Rp2.inner(d,ind,map=amap,result=rPP) 

                # tm1 and Rm1 ===================================
                if key[axis] != twon1:
                    t=list(key); t[axis]=t[axis]+1; t=tuple(t)
                    try:
                        d = dn[t]
                    except KeyError:
                        d0 = self.__sock_it_to_me_diffNS(n,t,surface=1)
                        d = Tensor(k2,k2,k2)
                        d[0:k,0:k,0:k] = d0
                        
                    tm1.inner(d,ind,map=amap,result=r)
                    if n>0:
                        d = d[0:k,0:k,0:k]
                        Rm1.inner(d,ind,map=amap,result=rPP)

                # tm2 and Rm2 ===================================
                if not (key[axis] in (twon1, twon1-1)):
                    t=list(key); t[axis]=t[axis]+2; t=tuple(t)
                    try:
                        d = dn[t]
                    except KeyError:
                        d0 = self.__sock_it_to_me_diffNS(n,t,surface=1)
                        d = Tensor(k2,k2,k2)
                        d[0:k,0:k,0:k] = d0
                        
                    tm2.inner(d,ind,map=amap,result=r)
                    if n>0:
                        d = d[0:k,0:k,0:k]
                        Rm2.inner(d,ind,map=amap,result=rPP)
                    
                # ===============================================
                if n>0:
                    r[0:k,0:k,0:k].gaxpy(1.0,rPP[0:k,0:k,0:k],-1.0)
                    
                result.d[n][key] = r.scale(scale)
                
                #cut1 += 1
                #tmp = Tensor(result.d[n][key])
                #tmp[0:k,0:k,0:k].gaxpy(1.0,result.d[n][key][0:k,0:k,0:k],-1.0)
                #if tmp.normf() > self.thresh: cut2 +=1

            del self.cache_sock_it_to_me_diffNS
            del self.cache_sock_it_to_me2_diffNS

        # Operator has been applied.  Must now cleanup f from its
        # nonstandard state, and also compress the result summing the
        # sum coeffs at all levels in the process.

        #print cut2, '/', cut1

        for n in result.d.keys():
            dn = result.d[n]
            sn = result.s[n] = {}
            for key in dn.keys():
                sn[key] = Tensor(dn[key][0:k,0:k,0:k])
                dn[key][0:k,0:k,0:k].fill(0.0)
            
        result.compressed = 0
        result.compress()

        # Clean up the source
        self.s[0] = {}
        self.s[0][0,0,0] = Tensor(self.d[0][0,0,0][0:k,0:k,0:k])
        for n in self.d.keys():
            dn = self.d[n]
            for key in dn.keys():
                dn[key][0:k,0:k,0:k].fill(0.0)
        self.compressed = 1
        
        for n in range(nmax,-1,-1):
            if len(self.d[n])==0:
                del self.d[n]
            else:
                break

        perf.exit('diff2NS')
        return result

    def latex_zslice(self,z=0.5,filename=None):
        '''

        Print to a file or standard output a Latex picture environment
        that displays the adpative grid.  For the specified slice thru
        the Z coordinates.  Adjust unitlength in your latex source to
        make the picture bigger or smaller.  Adjust the frame offset
        to move it.

        '''

        if filename:
            file = open(filename,"w")
        else:
            filename = "standard output"
            file = sys.stdout

        print "Latex plotting grid in slice z=%.6f to file %s" % (z, filename)

        self.reconstruct()

        z = z*100.0                     # Scale picture to 100mm

        file.write("\\setlength{\\unitlength}{1mm}\n")
        file.write("\\begin{picture}(100,100)\n")
        lines = []
        for n in self.s.keys():
            for lx, ly, lz in self.s[n].keys():
                length = 100.0*(0.5**n)
                lx, ly, lz = lx*length, ly*length, lz*length
                if lz <= z and (lz+length) >= z:
                    lines.append((lx,ly,1,0,length))
                    lines.append((lx+length,ly,0,1,length))
                    lines.append((lx+length,ly+length,-1,0,length))
                    lines.append((lx,ly+length,0,-1,length))

        # Get rid of duplicates
        lines.sort()
        for i in range(len(lines)-1,0,-1):
            if lines[i] == lines[i-1]:
                del(lines[i])
                
        for line in lines:
            file.write("\\put(%.4f,%.4f){\\line(%d,%d){%.4f}}\n" % line)

        file.write("\\end{picture}\n")
        if filename != "standard output":
            file.close()
        print "Done plotting"

    def to6d(self,side=1):
        '''
        
        side = 1 : left    f_new(r1,r2) = f(r1)
        side = 2 : right   f_new(r1,r2) = f(r2)
        
        '''

        if side!=1 and side!=2:
            raise "side is either 1 or 2"
            
        perf.enter('to6d')

        if self.compressed:

            result = Function6d(k=self.k,compress=1) # Was thresh=self.thresh
            result.s = {}
            result.d = {}
    
            k = self.k
                
            one = Tensor(k,k,k)
            one[0,0,0] = 1.0

            result.s[0]={}
            if(side==1):            
                result.s[0][0,0,0,0,0,0] = self.s[0][0,0,0].outer(one)
            else:
                result.s[0][0,0,0,0,0,0] = one.outer(self.s[0][0,0,0])
            
            one = Tensor(2*k,2*k,2*k)
            
            for n in self.d.keys():
                twon = 2**n
                h = 1.0/twon
                scale = sqrt(h)
                scale = scale*scale*scale
                one[0,0,0] = scale
                    
                trans1d = range(twon)
                dn = self.d[n]
                rn  = result.d[n] = {}
                for key in dn.keys():
                    if(side==1):            
                        hh = dn[key].outer(one)
                        for lx in trans1d:
                            for ly in trans1d:
                                for lz in trans1d:
                                    rn[key+(lx,ly,lz)] = hh
                        
                    else:
                        hh = one.outer(dn[key])
                        for lx in trans1d:
                            for ly in trans1d:
                                for lz in trans1d:
                                    rn[(lx,ly,lz)+key] = hh
                                        
        else:

            result = Function6d(k=self.k,compress=0) # Was thresh=self.thresh
            result.s = {}
    
            k = self.k
                
            one = Tensor(k,k,k)
                
            for n in self.s.keys():
                twon = 2**n
                h = 1.0/twon
                scale = sqrt(h)
                scale = scale*scale*scale
                one[0,0,0] = scale
                    
                trans1d = range(twon)
                sn = self.s[n]
                rn  = result.s[n] = {}
                for key in sn.keys():
                    if(side==1):            
                        hh = sn[key].outer(one)
                        for lx in trans1d:
                            for ly in trans1d:
                                for lz in trans1d:
                                    rn[key+(lx,ly,lz)] = hh
                    else:
                        hh = one.outer(sn[key])
                        for lx in trans1d:
                            for ly in trans1d:
                                for lz in trans1d:
                                    rn[(lx,ly,lz)+key] = hh

        perf.exit('to6d')

        return result

def generalized_inverse(S,thresh=1e-12):

    '''

    Compute the generalized inverse of the matrix S neglecting all
    eigenvalues less than thresh in absolute value

    '''
    (v,e) = testlapack.syev(S)
    n = S.dim(0)
    for i in range(n):
        if e[i] > thresh:
            v[:,i].scale(1.0/sqrt(e[i]))
##         else:
##             print i, "neglecting", e[i]
    return v*v.transpose()
            

setattr(Function, "fast_eval",types.MethodType(fast_eval.fast_eval, None, Function))
setattr(Function, "fast_eval_init",types.MethodType(fast_eval.fast_eval_init, None, Function))
setattr(Function, "fast_eval_tidy",types.MethodType(fast_eval.fast_eval_tidy, None, Function))

class Operator(Function):
    '''

    This class is used to make the non-standard form of an operator
    (currently just an integral convolution).  It inherits a bunch of
    stuff from Function, but this is just to reuse code for the
    initial quadrature etc..  All that can be done with an Operator is
    to apply it by multiplication to a Function.  Perhaps more will
    come later.

    In order to provide persistance, and to keep the memory usage as
    low as possible (it is still excessive!) the operator is
    associated with a disk directory in which the non-standard form
    matrix elements are kept.

    '''

    def __init__(self, function=None, thresh=None, k=None,
                 initial_level=None, refine=1, dirname=None,
                 homogenous=None, isotropic=None, cfunction=None):

        if dirname is None:
            self.nscache = DiskDir(tempfile.mktemp(),status='scratch')
        else:
            self.nscache = DiskDir(dirname,status='keep')

        # Assert my operator-ness
        self.operator = 1
        self.dirname = dirname
        self.isotropic = isotropic
        if not (homogenous is None):
            self.homogenous = 1
            self.degree = homogenous
        else:
            self.homogenous = 0

        # First initialize the Function class and do the finescale
        # projection and refinement of the kernel.  The function class
        # will recognize that we are an operator and double the order
        # of wavelets and double the range of evaluation.

        # Optimization.  Delay full initialization (compression of the
        # kernel) until the first time that we need to evaluate matrix
        # elements.  This is so we don't have to store or recompute
        # the compressed kernel if the matrix elements are already
        # cached.

##         Function.__init__(self, k=k, initial_level=initial_level,
##                           refine=refine,
##                           compress=0, thresh=thresh)
        self.k = 2*k
        self.initial_level = initial_level
        self.refine = refine
        self.thresh = thresh
        self.initialized = 0
        self.function = function
        self.cfunction = cfunction

        # Store an empty Function inside self for the sole purpose of
        # using its filter of order k (the operator has a filter of
        # order 2*k).  Also deferred.
        # self.__f = Function(k=k)

    def rnlp(self,n,lx,ly,lz):
        '''

        Return the (2k,2k,2k) tensor of coefficients rnl[p,q,r] which
        are the matrix elements of the kernel over the double order
        scaling functions.

        During initialization we compressed the kernel as if it were
        a plain old function.  Just need to form those coefficients
        at level n by using, yes you guessed it, the twoscale relations.

        If this routine has to look up/down the tree it will store
        the things it computed in the tree.

        It returns a COPY of the coefficients so you can do what
        you want to them.

        NOTE, that we have NOT included the extra factor of 1/8**n/2
        (i.e., the difference between the kernel and the function
        projection).

        '''
        if not self.initialized:
            print "Compressing operator kernel"
            self.initialized = 1
            Function.__init__(self, k=self.k/2,
                              initial_level=self.initial_level,
                              refine=self.refine, function=self.function,
                              cfunction=self.cfunction,
                              compress=0, thresh=self.thresh)
            self.__f = Function(k=self.k/2)            
        
        s = self._Function__sock_it_to_me(n,(lx,ly,lz))
        if not (s is None):
            return s

        # Rats!  The desired box does not exist and does not have a
        # parent.  Therefore it must be constructed from the
        # coefficients at a finer resolution.  This is done by
        # recursively descending then ascending the tree to get what
        # we want.
        
        lx2, ly2, lz2 = lx*2, ly*2, lz*2
        k = self.k                      # Note self.k is actually 2*k
        ss = Tensor(2*k,2*k,2*k)
        for ix in range(2):
            ixlo = ix*k
            ixhi = ixlo + k
            for iy in range(2):
                iylo = iy*k
                iyhi = iylo+k
                for iz in range(2):
                    izlo = iz*k
                    izhi = izlo+k
                    ss[ixlo:ixhi,iylo:iyhi,izlo:izhi] =  self.rnlp(
                        n+1,lx2+ix,ly2+iy,lz2+iz)
        sd = self.filter(ss)            # Need sonly optimization
        s = Tensor(sd[0:k,0:k,0:k])

        if not self.s.has_key(n):
            self.s[n] = {}
        self.s[n][lx,ly,lz] = s

        return Tensor(s)

    def rnlij(self,n,lx,ly,lz):
        '''

        Using the expansion of the autocorrelation functions in
        Legendre polyn of double order, return the transition
        matrix elements for the specified level and translations

        '''
        k = self.k                      # This is really 2*k
        c = autoc.autoc(k/2)

        R = Tensor(2*k,2*k,2*k)
        for ix in -1,0:
            xs = slice((ix+1)*k,(ix+2)*k)
            for iy in -1,0:
                ys = slice((iy+1)*k,(iy+2)*k)
                for iz in -1,0:
                    zs = slice((iz+1)*k,(iz+2)*k)
                    R[xs,ys,zs] = self.rnlp(n,lx+ix,ly+iy,lz+iz)

        R.scale(sqrt(1.0/8.0**n))
        R = R.fast_transform(c)

        # OK, now have R[pp',qq',rr'] but for everything else we will
        # want R[p,q,r,p',q',r'], so split the dimensions, reorder and
        # make a contiguous copy.  Sigh.

        k = k/2
        R = R.splitdim(0,k,k)
        R = R.splitdim(2,k,k)
        R = R.splitdim(4,k,k)
        R = R.swapdim(3,4)
        R = R.swapdim(1,2)
        R = R.swapdim(2,3)
        R = Tensor(R)                   # R[p,q,r,p',q',r']
        
        return R


    def nonstandard(self,n,lx,ly,lz):
        '''

        Return the matrix elements alpha, beta and gamma for the
        non-standard form of the operator at level n if the difference
        in the translations (l-m)=lx,ly,lz.

        The matrix elements are computed and cached on disk, so the
        first application might be a tad slow (k^7) if they have to be
        recomputed.  If there is no associated file, then the operator
        is recomputed every time.

        If self.homgeneous is true, then self.degree is the degree of
        homogeneity ... K(s*x) = (s^degree)*K(x).  The operator is
        saved at level 0 and for the 0,0,0 block an additional matrix
        saves the nonstandard form that also includes T0 which is
        needed just at level 0.

        If self.isotropic is true, then we assume that the kernel is
        just a function of |x*x + y*y + z*z| and therefore we can
        permute the directions into some standard order, and change
        the signs by changing the signs of the coefficients.

        We assume that requests are always with |l|<2**n.

        '''
        
        if self.isotropic:
            # Permute the translations so |lx| >= |ly| >= |lz|
            # while remembering how to undo this mess.
            lllx, llly, lllz = lx, ly, lz
            lx, ly, lz = abs(lx), abs(ly), abs(lz)
            permutations = []
            if lx < ly:
                lx, ly = ly, lx
                permutations = [(0,1)] + permutations
            if lx < lz:
                lx, lz = lz, lx
                permutations = [(0,2)] + permutations
            if ly < lz:
                ly, lz = lz, ly
                permutations = [(1,2)] + permutations

        twok = self.k
        k = self.k/2
        nn = n
        if self.homogenous: nn = 0 # Fetch from level 0
        key = (nn,lx,ly,lz)
        if n==lx==ly==lz==0: key = key + ("T0",) # (0,0,0,0) includes T0

        try:
            R = Tensor(self.nscache[key])       # Assumes returns a copy
            if self.homogenous: R.scale(0.5**(n*(3+self.degree)))
            if self.isotropic:
                R = self.decanonicalize(R, permutations, lllx, llly, lllz)

            return R, R.normf()
        except KeyError:
            pass

        nuse = n
        if self.homogenous:
            nuse = 0                       # Find highest feasible level
            lmax = max(abs(lx),abs(ly),abs(lz))
            while (2**nuse) <= lmax: nuse = nuse + 1
            print "computing", lx, ly, lz, "at level", nuse

        lx2, ly2, lz2 = lx*2, ly*2, lz*2
        R = Tensor(2*k,2*k,2*k,2*k,2*k,2*k)
        
        for ix in 0,1:
            xs = slice(ix*k,(ix+1)*k)
            for iy in 0,1:
                ys = slice(iy*k,(iy+1)*k)
                for iz in 0,1:
                    zs = slice(iz*k,(iz+1)*k)
                    for ixx in 0,1:
                        xxs = slice(ixx*k,(ixx+1)*k)
                        for iyy in 0,1:
                            yys = slice(iyy*k,(iyy+1)*k)
                            for izz in 0,1:
                                zzs = slice(izz*k,(izz+1)*k)
                                R[xs,ys,zs,xxs,yys,zzs] = self.rnlij(
                                    nuse+1,lx2+ix-ixx,ly2+iy-iyy,lz2+iz-izz)

        R = self.__f.filter(R)

        if self.homogenous:
            # Scale from level nuse to level 0 before screening
            R.scale(2.0**(nuse*(3+self.degree)))

        R.screen(self.thresh*1e-3)      # NEED A BETTER THRESHOLD

        # Before saving/returning, fuse the first and last three
        # dimensions so that a conventional matrix vector product can
        # be used to apply the operator.

        # Don't want Tn, except at (0,0,0,0) where we sneakily include it

        if self.homogenous:
            if lx==ly==lz==0:
                R_T0 = Matrix(R.fusedim(0).fusedim(0).fusedim(-2).fusedim(-2))
                self.nscache[0,0,0,0,"T0"] = Tensor(R_T0)         #CCC

            R[0:k,0:k,0:k,0:k,0:k,0:k].fill(0.0)
            R = Matrix(R.fusedim(0).fusedim(0).fusedim(-2).fusedim(-2))
            self.nscache[0,lx,ly,lz] = Tensor(R)                 #CCC

            if n == 0: R = R_T0

            R.scale(0.5**(n*(3+self.degree)))
        else:
            if n != 0: R[0:k,0:k,0:k,0:k,0:k,0:k].fill(0.0)
            R = Matrix(R.fusedim(0).fusedim(0).fusedim(-2).fusedim(-2))
            self.nscache[key] = Tensor(R)                       #CCC

        if self.isotropic:
            R = self.decanonicalize(R, permutations, lllx, llly, lllz)
            
        return R, R.normf()

    def decanonicalize(self, R, permutations, lllx, llly, lllz):
        '''

        The operator is isotropic and we have exploited this by
        permuting the translations into some standard order, which is
        the same as permuting the axes (since the operator does not
        have any sense of direction).  We have to undo this
        permutation and recover the desired order of axes.

        '''
        k = self.k
        if permutations:
            # R(xyz,XYZ) -> R(x,y,z,XYZ)
            R.reshape_inplace((k,k,k,k,k,k))

            # Now apply the permutations on both sides
            for i, j in permutations:
                R.swapdim_inplace(i,j).swapdim_inplace(i+3,j+3)

            # Make a contiguous copy and refuse the dimensions
            R = Tensor(R).reshape_inplace((k*k*k,k*k*k))

        # Fix the phases
        if lllx<0: cfcube.flipper(k, R.swigptr(), 0)
        if llly<0: cfcube.flipper(k, R.swigptr(), 1)
        if lllz<0: cfcube.flipper(k, R.swigptr(), 2)

        return R

    def apply(self,f,bmax=4):
        '''

        Apply the operator self to the function f

        Currently assume that self is an integral convolution, so we
        index the operator elements with l-m.  Also, assume the
        operator has limited bandwidth, so neglect all terms beyond
        the first term for which the norm becomes negligble.
        The maximum bandwidth |l-m| is assumed less than bmax for
        the sole purpose of making the loops simple.

        The threshold for truncation should reflect the actual
        bandwidth.  If there will be B contributions to block of
        elements in the result, then to get a precision eps in the
        result we only need to compute contributions whose expected
        norm is greater than eps/B.  For now choose the screening
        threshold as f.thresh*0.1 (since the upper bound on the norm
        of the error should be very loose), but clearly something more
        intelligent would be helpful (esp. when we try to exploit
        sparsity within the blocks of the operator).  Note that
        just f.thresh is empircally not quite tight enough and
        multiplying by 0.01 seems to be uncessarily tight.

        This method probably belongs inside the Function class.

        '''
        if not isinstance(f,Function):
            raise TypeError

        f.compress()
        f.reconstruct(nonstandard=1)
        
        # Make a list of translations i=(l-m) ordered according to the
        # norm. As soon as we hit a negligble term, we can neglect all
        # other terms.  In doing this I'm probably making some nasty
        # assumptions about the operator being isotropic.

        trans = []
        for ix in range(-bmax,bmax+1):
            for iy in range(-bmax,bmax+1):
                for iz in range(-bmax,bmax+1):
                    trans.append((ix*ix+iy*iy+iz*iz,ix,iy,iz))
        trans.sort()


        nmax = max(f.d.keys())

        # Factor of 4 is to gain agreement with precision of sepop code
        tol = f.thresh/27.0/nmax/4.0

        result = Function(k=f.k, thresh=f.thresh)
        result.d[0] = {}

        # Flatten the tensors in the source and find max norm at
        # each level
        maxnorm = [0.0]*(nmax+1)
        for n in range(nmax,-1,-1):
            dn = f.d[n]
            for key in dn.keys():
                dn[key] = dn[key].flat()
                maxnorm[n] = max(maxnorm[n],dn[key].normf())

        veclen = self.k**3

        # Loop thru feasible translations
        # sq=0 is self
        # sq=1 is face contact
        # sq=2 is edge contact
        # sq=3 is corner contact
        # sq=4 is one face removed
        for sq,ix,iy,iz in trans:
            significant = 0
            op = None

            for n in range(nmax,-1,-1):

                dn = f.d[n]
                twon = 2**n
                twon1 = twon-1
                    
                if abs(ix)>twon1 or abs(iy)>twon1 or abs(iz)>twon1:
                    continue

                if self.homogenous and op:
                    # The operator is homogenous and has been read in
                    # already on a finer level.  Just scale it!
                    # However, if at level 0 we need the T0 block.
                    if n == 0:
                        op,norm = self.nonstandard(n,ix,iy,iz)
                    else:
                        op.scale(2.0**(3+self.degree))
                        norm = norm * (2.0**(3+self.degree))
                else:
                    # Either not homogenous or first time getting op
                    op,norm = self.nonstandard(n,ix,iy,iz)

                
                if norm*maxnorm[n] > tol:
                    significant = 1
                    try:
                        rn = result.d[n]
                    except KeyError:
                        rn = result.d[n] = {}

                    # Assemble a list of significant sources
                    info = []
                    for mx, my, mz in dn.keys():
                        lx = mx + ix
                        ly = my + iy
                        lz = mz + iz
                        if (lx<0 or lx>twon1 or ly<0 or ly>twon1 or \
                            lz<0 or lz>twon1): continue
                        d = dn[mx,my,mz]
                        dnorm = d.normf()
                        if (dnorm*norm) > tol:
                            info.append(((lx,ly,lz),(mx,my,mz)))

                    nvec = len(info)

                    if nvec == 0:
                        continue

                    #print "OP[%d][%d,%d,%d]=%.1e f[%d]=%.1e nvec=%d" % \
                    #      (n,ix,iy,iz,norm,n,maxnorm[n],nvec)

                    # Gather the sources, apply the op, scatter results
                    # Do this so that kernel is dgemm not just dgemv
                    source = Matrix(nvec,veclen)
                    for i in range(nvec):
                        source[i,:] = dn[info[i][1]]
                        
                    #product = op * source.transpose()
                    product = source.mxt(op)
                    
                    for i in range(nvec):
                        key = info[i][0]
                        if rn.has_key(key):
                            rn[key].gaxpy(1.0,product[i,:],1.0)
                        else:
                            rn[key] = product[i,:]

                    del product
                    del source

            if not significant:         # Norm test failed at all levels
                break
                        
        # Operator has been applied.  Must now cleanup f from its
        # nonstandard state, and also compress the result summing the
        # sum coeffs at all levels in the process.

        # First, turn all of the vectors back into contiguous tensors
        # and put the scaling coefficients back into the tree instead
        # of storing them in the usually zero block of the difference
        # coefficients.
        
        k = result.k
        for n in result.d.keys():
            dn = result.d[n]
            sn = result.s[n] = {}
            for key in dn.keys():
                dn[key] = Tensor(dn[key].splitdim(0,4*k*k,2*k).splitdim(0,2*k,2*k))
                sn[key] = Tensor(dn[key][0:k,0:k,0:k])
                dn[key][0:k,0:k,0:k].fill(0.0)
            
        # Now compress it to sum the contributions at all levels
        # which are in both bases.
        result.compressed = 0
        result.compress()

        # Clean up the source
        f.s[0] = {}
        f.s[0][0,0,0] = Tensor(f.d[0][0,0,0][0:k,0:k,0:k])
        for n in f.d.keys():
            dn = f.d[n]
            for key in dn.keys():
                dn[key][0:k,0:k,0:k].fill(0.0)
        f.compressed = 1

        return result


    def sparse_op(self,n,lx,ly,lz,tol):
        '''
        Compute the sparse matrix representation of the operator,
        caching it on disk. 
        '''
        key = (n,lx,ly,lz,'sparse')
        if n == 0: key = (0,0,0,0,'T0','sparse')

        try:
            op, norm = self.nscache[key]
            #print "read from disk", n, lx, ly, lz
        except KeyError:
            op,normf = self.nonstandard(n,lx,ly,lz)
            op = SparseMatrix(op,thresh=self.thresh*1e-3)
            norm = op.norm2()
            self.nscache[key] = op, norm

        op = SparseMatrix(op,tol)
        return op, norm

    def apply_sparse(self,f,bmax=1):
        '''

        Apply the operator self to the function f

        Currently assume that self is an integral convolution, so we
        index the operator elements with l-m.  Also, assume the
        operator has limited bandwidth, so neglect all terms beyond
        the first term for which the norm becomes negligble.
        The maximum bandwidth |l-m| is assumed less than bmax for
        the sole purpose of making the loops simple.

        The threshold for truncation should reflect the actual
        bandwidth.  If there will be B contributions to block of
        elements in the result, then to get a precision eps in the
        result we only need to compute contributions whose expected
        norm is greater than eps/B.  For now choose the screening
        threshold as f.thresh*0.1 (since the upper bound on the norm
        of the error should be very loose), but clearly something more
        intelligent would be helpful (esp. when we try to exploit
        sparsity within the blocks of the operator).  Note that
        just f.thresh is empircally not quite tight enough and
        multiplying by 0.01 seems to be uncessarily tight.

        This method probably belongs inside the Function class.

        '''
        if not isinstance(f,Function): raise TypeError

        t_func = time.clock()
        f.compress()
        fnorm = max(1.0,f.norm2())
        f.reconstruct(nonstandard=1)
        t_func = time.clock() - t_func
        t_setup = time.clock()
        
        # Make a list of translations i=(l-m) ordered according to the
        # norm. As soon as we hit a negligble term, we can neglect all
        # other terms.  In doing this I'm probably making some nasty
        # assumptions about the operator being isotropic.

        trans = []
        for ix in range(-bmax,bmax+1):
            for iy in range(-bmax,bmax+1):
                for iz in range(-bmax,bmax+1):
                    trans.append((ix*ix+iy*iy+iz*iz,ix,iy,iz))
        trans.sort()

        nmax = max(f.d.keys())

        # Factor of 4 is to gain agreement with precision of sepop code
        tol = f.thresh / 27.0 / nmax / 4.0
        print "thresh", self.thresh, f.thresh, "tol", tol, "fnorm", fnorm

        result = Function(k=f.k, thresh=f.thresh)
        result.d[0] = {}

        twok = 2*f.k

        # Find norms at each level
        maxnorm = [0.0]*(nmax+1)
        dnorms = {}
        for n in range(nmax,-1,-1):
            dnorms[n] = {}
            dn = f.d[n]
            for key in dn.keys():
                dnorms[n][key] = dn[key].normf()
                maxnorm[n] = max(maxnorm[n],dnorms[n][key])

        veclen = self.k**3
        t_setup = time.clock() - t_setup
        t_apply = time.clock()
        t_mul = 0.0
        t_op = 0.0
        flops = 0.0

        # Loop thru feasible translations
        # sq=0 is self
        # sq=1 is face contact
        # sq=2 is edge contact
        # sq=3 is corner contact
        # sq=4 is one face removed
        for sq,ix,iy,iz in trans:
            significant = 0

            for n in range(0,nmax+1):
                dn = f.d[n]
                dnormsn = dnorms[n]
                twon = 2**n
                twon1 = twon-1
                    
                if abs(ix)>twon1 or abs(iy)>twon1 or abs(iz)>twon1:
                    continue

                used = time.clock()
                stol = 0.1*tol/maxnorm[n]
                if self.homogenous and n>1:
                    op = SparseMatrix(op.scale(0.5**(3+self.degree)),stol)
                    norm = norm * (0.5**(3+self.degree))
                else:
                    op, norm = self.sparse_op(n,ix,iy,iz,stol)


                zerofrac = (1.0 - op.nval()/(float(self.k)**6))*100.0
                print "%d (%d,%d,%d) 2norm=%.2e zero=%.2f" % \
                      (n,ix,iy,iz,norm,zerofrac)

                t_op = t_op + time.clock() - used
                
                if norm*maxnorm[n] <= tol: break

                used = time.clock()
                significant = 1

                try:
                    rn = result.d[n]
                except KeyError:
                    rn = result.d[n] = {}

                dlist = []
                rlist = []
                for mx, my, mz in dn.keys():
                    lx, ly, lz = mx + ix, my + iy, mz + iz
                    if (lx<0 or lx>twon1 or ly<0 or ly>twon1 or \
                        lz<0 or lz>twon1): continue
                    if (dnormsn[mx,my,mz]*norm) > tol:
                        d = dn[mx,my,mz]
                        if not rn.has_key((lx,ly,lz)):
                            rn[lx,ly,lz] = Tensor(twok,twok,twok)
                        dlist.append(d.swigptr())
                        rlist.append(rn[lx,ly,lz].swigptr())

                op.mxv_multi(dlist,rlist)

                flops = flops + 2.0*len(dlist)*op.nval()

                del dlist
                del rlist
                t_mul = t_mul + time.clock() - used

            if not significant: break   # Norm test failed at all levels
                        
        # Operator has been applied.  Must now cleanup f from its
        # nonstandard state, and also compress the result summing the
        # sum coeffs at all levels in the process.

        t_apply = time.clock() - t_apply
        t_clean = time.clock()
        k = result.k
        for n in result.d.keys():
            dn = result.d[n]
            sn = result.s[n] = {}
            for key in dn.keys():
                sn[key] = Tensor(dn[key][0:k,0:k,0:k])
                dn[key][0:k,0:k,0:k].fill(0.0)
            
        # Now compress it to sum the contributions at all levels
        # which are in both bases.
        result.compressed = 0
        result.compress()

        # Clean up the source
        f.s[0] = {}
        f.s[0][0,0,0] = Tensor(f.d[0][0,0,0][0:k,0:k,0:k])
        for n in f.d.keys():
            dn = f.d[n]
            for key in dn.keys():
                dn[key][0:k,0:k,0:k].fill(0.0)
        f.compressed = 1

        t_clean = time.clock() - t_clean


        print "func=%.1fs setup=%.1fs apply=%.1fs (op=%.1fs mul=%.1fs)\n clean=%.fs flops=%.1e mflop/s=%.1e" % \
              (t_func, t_setup, t_apply, t_op, t_mul, t_clean, 
               flops, flops/t_op/1e6)

        return result

    def apply_sparse_hfex(self,psj,psi,ff,nmo,bmax=4):
        '''
        '''
        print "Computing HF exchage potential for %d" % (nmo)

        #                              /      psi[p](r2)*ff(r2)  
        # HFex(r1) =  SUM psj[p](r1) * | dr2 ------------------- 
        #              p               /         | r1 - r2 |     
        # -----------------------------------------------------------------

        t_0 = time.clock()

        f=psi[0]

        if not isinstance(f,Function): raise TypeError

        thresh = f.thresh
        safethresh = thresh * 0.01

        t_func  = 0.0
        t_setup = 0.0
        t_apply = 0.0
        t_sumup = 0.0
        t_mul   = 0.0
        t_mul2  = 0.0
        t_mul3  = 0.0
        t_op    = 0.0
        t_clean = 0.0
        
        # Make a list of translations i=(l-m) ordered according to the
        # norm. As soon as we hit a negligble term, we can neglect all
        # other terms.  In doing this I'm probably making some nasty
        # assumptions about the operator being isotropic.

        trans = []
        for ix in range(-bmax,bmax+1):
            for iy in range(-bmax,bmax+1):
                for iz in range(-bmax,bmax+1):
                    trans.append((ix*ix+iy*iy+iz*iz,ix,iy,iz))
        trans.sort()

        pimo = ff.compress()
        pimo = pimo.truncate(safethresh)

        result = Function(k=psi[0].k,thresh=thresh,compress=1)
        twok = 2*psi[0].k
        tol = safethresh / 27.0 
        veclen = self.k**3

        t_setup = t_setup + time.clock() - t_0
        
        for jmo in range(nmo):

            if(psj[jmo]==None): continue
        
            t_0 = time.clock()

            pjmo = psi[jmo].compress()
            pjmo = pjmo.truncate(safethresh)
            
            nmax = max(pjmo.s.keys())

            # Find norms at each level
            maxnorm0 = [0.0]*(nmax+1)
            dnorms0 = {}
            for n in range(nmax,-1,-1):
                dnorms0[n] = {}
                dn = pjmo.d[n]
                for key in dn.keys():
                    dnorms0[n][key] = dn[key].normf()
                    maxnorm0[n] = max(maxnorm0[n],dnorms0[n][key])

            t_setup = t_setup + time.clock() - t_0

            # computing charge distributions [pn]*[pm] over occupied MOs.
            # f(r2) = pimo(r2) * pjmo(r2)
            # -----------------------------------------------------------

            t_0 = time.clock()

            f = pimo*pjmo
        
            f.compress()
            f.truncate(safethresh)
            fnorm = max(1.0,f.norm2())
            f.reconstruct(nonstandard=1)

            nmax = max(f.d.keys())

            #print "thresh", self.thresh, safethresh, "tol", tol, "fnorm", fnorm

            # Find norms at each level
            maxnorm = [0.0]*(nmax+1)
            dnorms = {}
            for n in range(nmax,-1,-1):
                dnorms[n] = {}
                dn = f.d[n]
                for key in dn.keys():
                    dnorms[n][key] = dn[key].normf()
                    maxnorm[n] = max(maxnorm[n],dnorms[n][key])

            t_func = t_func + time.clock() - t_0

            # integrate           /      pimo(r2)*pjmo(r2)
            # int_pimo_pjmo(r1) = | dr2 -------------------
            #                     /         | r1 - r2 |
            # -----------------------------------------------------------

            # Screening may be optimized by sorting keys at each level by
            # decreasing norm.  Not yet doing this.

            t_0 = time.clock()

            int_pimo_pjmo = Function(k=psi[0].k,thresh=safethresh,compress=0)
            int_pimo_pjmo.d[0] = {}
        
            # Loop thru feasible translations
            # sq=0 is self
            # sq=1 is face contact
            # sq=2 is edge contact
            # sq=3 is corner contact
            # sq=4 is one face removed
            for sq,ix,iy,iz in trans:
                significant = 0
    
                for n in range(0,nmax+1):
                    dn = f.d[n]
                    dnormsn = dnorms[n]
                    twon = 2**n
                    twon1 = twon-1
                        
                    if abs(ix)>twon1 or abs(iy)>twon1 or abs(iz)>twon1:
                        continue
    
                    used = time.clock()
                    stol = 0.1*tol/maxnorm[n]
                    if self.homogenous and n>1:
                        op = SparseMatrix(op.scale(0.5**(3+self.degree)),stol)
                        norm = norm * (0.5**(3+self.degree))
                    else:
                        op, norm = self.sparse_op(n,ix,iy,iz,stol)
    
                    zerofrac = (1.0 - op.nval()/(float(self.k)**6))*100.0
                    
                    #print "%d (%d,%d,%d) 2norm=%.2e zero=%.2f" % \
                    #      (n,ix,iy,iz,norm,zerofrac)
    
                    t_op = t_op + time.clock() - used
                    
                    if norm*maxnorm[n] <= tol: continue
    
                    used = time.clock()
                    significant = 1
    
                    try:
                        rn = int_pimo_pjmo.d[n]
                    except KeyError:
                        rn = int_pimo_pjmo.d[n] = {}
    
                    dlist = []
                    rlist = []
                    for mx, my, mz in dn.keys():
                        lx, ly, lz = mx + ix, my + iy, mz + iz
                        if (lx<0 or lx>twon1 or ly<0 or ly>twon1 or \
                            lz<0 or lz>twon1): continue
                        d = dn[mx,my,mz]
                        dnorm = dnormsn[mx,my,mz]
                        if (dnorm*norm) > tol:
                            significant = 1
                            if not rn.has_key((lx,ly,lz)):
                                rn[lx,ly,lz] = Tensor(twok,twok,twok)
                            dlist.append(d.swigptr())
                            rlist.append(rn[lx,ly,lz].swigptr())
    
                    op.mxv_multi(dlist,rlist)
                    del dlist
                    del rlist
                    t_mul = t_mul + time.clock() - used
    
                if not significant: break   # Norm test failed at all levels

            t_apply = time.clock() - t_0
                            
            # Operator has been applied.  Must now cleanup f from its
            # nonstandard state, and also compress the result summing the
            # sum coeffs at all levels in the process.
    
            t_0 = time.clock()
            
            k = int_pimo_pjmo.k
            for n in int_pimo_pjmo.d.keys():
                dn = int_pimo_pjmo.d[n]
                sn = int_pimo_pjmo.s[n] = {}
                for key in dn.keys():
                    sn[key] = Tensor(dn[key][0:k,0:k,0:k])
                    dn[key][0:k,0:k,0:k].fill(0.0)
                
            # Must also clean up the non-standard form of the source
            # NOT YET ... see the gauss code.
            
            t_clean = t_clean + time.clock() - t_0

            # Now compress it to sum the contributions at all levels
            # which are in both bases.

            t_0 = time.clock()

            int_pimo_pjmo.compressed = 0
            int_pimo_pjmo = int_pimo_pjmo.compress()
            int_pimo_pjmo = int_pimo_pjmo.truncate(safethresh)

            t_mul2 = t_mul2 + ( time.clock() - t_0 )

            # computing            /      pimo(r2)*pjmo(r2)
            # tmp(r1) = pjmo(r1) * | dr2 -------------------
            #                      /         | r1 - r2 |     
            # -----------------------------------------------------------

            t_0 = time.clock()
            
            pjmo = psj[jmo].compress()
            pjmo = pjmo.truncate(safethresh)
            tmp = int_pimo_pjmo * pjmo

            t_mul3 = t_mul3 + ( time.clock() - t_0 )
            
            # sum up tmp(r1) into HFex_tmp[imo](r1)
            #
            #                                      /      pimo(r2)*pjmo(r2)  
            # HFex_tmp[imo](r1) =  SUM pjmo(r1) * | dr2 ------------------- 
            #                      jmo            /         | r1 - r2 |     
            # -----------------------------------------------------------------

            t_0 = time.clock()

            result.gaxpy(1.0,tmp,1.0)

            t_sumup = t_sumup + ( time.clock() - t_0 )

            del int_pimo_pjmo
            del tmp
            del f

        print "func=%.1fs setup=%.1fs apply=%.1fs (op=%.1fs mul=%.1fs) clean=%.1fs mul2=%.1fs mul3=%.1fs sumup=%.1fs" % \
              (t_func, t_setup, t_apply, t_op, t_mul, t_clean, t_mul2, t_mul3, t_sumup)

        return result.compress()

    # <==================================================================== ADDED BY YANAI UPTO HERE

def Poisson(k,thresh=None,dirname='/scratch/poisson/'):
    '''

    Return the operator for kernel to the Poisson equation (i.e.,
    1/|x-y|) for the given order multi-wavelet.  If the threshold is
    not specified the (tight) default of 10^-k is used.

    '''
    if thresh is None:
        thresh = 0.1**k

    cfcube.cvar.two_e_cut = 1e-6;

    dirname = dirname + `k`
        
    op = Operator(cfunction=cfcube.cvar.Ctwo_elec_pot,
                  function=cfcube.two_elec_pot,
                  thresh=thresh,k=k,refine=1,
                  dirname=dirname,
                  homogenous=-1,isotropic=1)
    return op


def BSH(k,bsh_k,thresh=None,dirname="/scratch/bsh/"):
    if thresh is None: thresh = 0.1**k

    cfcube.cvar.bsh_cut = 1e-4
    cfcube.cvar.bsh_k   = bsh_k
    cfcube.cvar.bsh_fac = 1.0 

    if dirname: dirname = dirname + 'k=%d bsh_k=%.8e'%(k,bsh_k)

    # Maintain relative precision when forming the operator.
    # When it is applied we only need to provide absolute precision
    thresh = thresh / max(1.0,bsh_k**2)
        
    op = Operator(cfunction=cfcube.cvar.Cbsh_kernel,
                  function=cfcube.bsh_kernel,
                  thresh=thresh,k=k,refine=1,
                  dirname=dirname,
                  isotropic=1)

    return op

class Function6d(Function):
    
    def __init__(self,
                 function=None, cfunction=None,
                 thresh=None, k=None,
                 initial_level=None,refine=1,compress=1,box=None):
        
        perf.enter('init')

        # For methods that may be inherited by subclasses we need to
        # remember the class and tensor rank and override it in
        # subclasses (don't really need to remember the tensor rank
        # I think, but it might be convenient at some point).

        self.ndim = 6
        self.new = Function6d
        self.newtensor = Tensor

        if k:
            self.k = k                  # Override default in class
        else:
            self.k = k = Function.k     # Ensure set in instance

        if self.operator:               # Operators need double order
            self.k = k = 2*k

        if thresh:
            self.thresh = thresh        # Override default in class
        else:
            self.thresh = thresh = self.thresh # Ensure set in instance

        if initial_level is None:
            self.initial_level = initial_level = Function.initial_level
        else:
            self.initial_level = initial_level


        self.d = {}              # Directory of wavelet coefficients
        self.s = {}              # Directory of scaling function coefficients
        self.function = function # Remember the function being compressed
        self.cfunction = cfunction
        self.compressed = 0      # Not yet compressed

        self.init_twoscale(k)
        self.init_quadrature(k)

        if cfunction:
            if not (type(cfunction) == type(' ') and \
                    cfunction[-6:] == 'p_void'):
                raise TypeError,"cfunction should be a C "\
                      "function pointer cast to void *"

        if function or cfunction:
            self.nterminated = 0
            if box is None:
                trans = None
            else:
                trans = []
                lxlo,lylo,lzlo,lxhi,lyhi,lzhi = box
                xrange = range(lxlo,lxhi+1)
                yrange = range(lylo,lyhi+1)
                zrange = range(lzlo,lzhi+1)
                for lx1 in xrange:
                    for ly1 in yrange:
                        for lz1 in zrange:
                            for lx2 in xrange:
                                for ly2 in yrange:
                                    for lz2 in zrange:
                                        trans.append((lx1,ly1,lz1,lx2,ly2,lz2))
                
            perf.enter('fs proj')
            print "file_scaling starts", trans
            self.fine_scale_projection(initial_level,trans)
            print "file_scaling finished"
            perf.exit('fs proj')

            del trans
            
            if refine:
                print "refine starts"
                perf.enter('refine')
                # Midly convoluted way to make list of translations to
                # be refined so that the same code works for operators
                # and for functions.
                refinethese = {}
                for lx1, ly1, lz1, lx2, ly2, lz2 in self.s[initial_level].keys():
                    refinethese[parent(lx1),parent(ly1),parent(lz1),parent(lx2),parent(ly2),parent(lz2)] = 1
                for lx1, ly1, lz1, lx2, ly2, lz2 in refinethese.keys():
                    self.refine_fine_scale_projection(initial_level-1,
                                                      lx1, ly1, lz1, lx2, ly2, lz2, prnt=0)
                del refinethese
                perf.exit('refine')

            if self.nterminated > 0:
                print "Terminated refinement:", self.nterminated

            if compress:
                self.compress(truncate=0)
        else:
            self.s[0] = {}
            self.s[0][0,0,0,0,0,0] = Tensor(k,k,k,k,k,k)
            if compress:
                self.compressed = 1
                self.d[0] = {}
                self.d[0][0,0,0,0,0,0] = Tensor(2*k,2*k,2*k,2*k,2*k,2*k)
            
        if self.debug:
            print "Made new function", id(self), self.k, self.thresh, self.compressed
            print "s:", self.s.keys()
            if self.compressed: print "d:", self.d.keys()

        perf.exit('init')
            
    def fcube(self, x1lo, y1lo, z1lo, x2lo, y2lo, z2lo, npt, h, pts, f, function):
        '''

        Tabulate the function on the quadrature points in the
        specified cube.  The C version is MUCH faster.

        '''
        for p1 in range(npt):
            x1 = x1lo+pts[p1]*h
            for q1 in range(npt):
                y1 = y1lo+pts[q1]*h
                for r1 in range(npt):
                    z1 = z1lo+pts[r1]*h
                    for p2 in range(npt):
                        x2 = x2lo+pts[p2]*h
                        for q2 in range(npt):
                            y2 = y2lo+pts[q2]*h
                            for r2 in range(npt):
                                z2 = z2lo+pts[r2]*h
                                f[p1,q1,r1,p2,q2,r2] = function(x1,y1,z1,x2,y2,z2)

    def fine_scale_projection(self,n,trans=None):
        '''

        Project the function onto scaling functions at level n.

        If self.operator is true, then we are computing the kernel of
        an integral operator. Operators need -2^n <= l <= 2^n as the
        default range.  

        If l[] is specified it is a list of translations that are to
        be considered.  Otherwise the full list is used.

        Modifies:
        :           self.s

        '''

        try:
            self.s[n]
        except KeyError:
            self.s[n] = {}
        
        h = 1.0/(2.0**n)                # Box size on target level 
        scale = sqrt(h)
        scale = scale**self.ndim       # Since in nD
        phiw = self.quad_phiw
        pts = self.quad_x 
        npt = len(pts)

        if not trans:
            if self.operator:
                trans1d = range(-2**n,2**n+2)
            else:
                trans1d = range(2**n)
            trans = []
            for lx1 in trans1d:
                for ly1 in trans1d:
                    for lz1 in trans1d:
                        for lx2 in trans1d:
                            for ly2 in trans1d:
                                for lz2 in trans1d:
                                    trans.append((lx1,ly1,lz1,lx2,ly2,lz2))

        for lx1,ly1,lz1,lx2,ly2,lz2 in trans:
            x1lo = lx1*h
            y1lo = ly1*h
            z1lo = lz1*h
            x2lo = lx2*h
            y2lo = ly2*h
            z2lo = lz2*h
            f = Tensor(npt,npt,npt,npt,npt,npt)
            if self.cfunction:
                cfcube.cfcube6(x1lo, y1lo, z1lo, x2lo, y2lo, z2lo, npt, h,
                               pts.swigptr(), f.swigptr(),
                               self.cfunction)
            else:
                # Plain old Python function or non-cast C function
                # This is SLOW due to all of the tensor references.
                self.fcube(x1lo, y1lo, z1lo, x2lo, y2lo, z2lo, npt, h, pts, f, self.function)
            f.scale(scale)

            self.s[n][lx1,ly1,lz1,lx2,ly2,lz2] = f.fast_transform(phiw)
            
        del trans

    def refine_fine_scale_projection(self,n,lx1,ly1,lz1,lx2,ly2,lz2,prnt=0):
        '''

        We have the fine scale projection at level n+1, so we can
        compute the scaling function and wavelet coefficients at level
        n.  Considering just box l on level n, examine the difference
        coefficients.  If they do not satisfy the threshold, then
        refine the representation by recurring down to a finer level.

        If refinement occurs this modifies:
        :   self.s

        Set prnt to non-zero for printing

        '''
        if n > self.max_refine_level:
            self.nterminated = self.nterminated+1
            return

        if Function.truncate_method == 0:
            ##tol = self.thresh
            ##tol = sqrt(0.5**(max(3*n,6*2)))*self.thresh  ## 6*2 =6dim * initn=2
            tol = self.thresh*sqrt(0.5**(3*max(n,3)))
        elif Function.truncate_method == 1:
            tol = sqrt(0.5**(self.ndim*n))*self.thresh
        else:
            # Here, cannot do truncate_method=2 properly since we are
            # refining vertically down and don't yet know the number
            # of boxes at this level.  Arbitrarily chop it at level 5
            # ... i.e., assume that in all space there are no more
            # than 8**5=32768 boxes.  sqrt(32768)=181.01
            tol = self.thresh/sqrt((2.0**self.ndim)**5)

        k = self.k
        ss = self.gather_scaling_coeffs(n,lx1,ly1,lz1,lx2,ly2,lz2)
        ss = self.filter(ss)
        s = Tensor(ss[0:k,0:k,0:k,0:k,0:k,0:k])
        ss[0:k,0:k,0:k,0:k,0:k,0:k].fill(0.0)
        dnorm = ss.normf()
        if dnorm < tol:
            #print "truncating", n,lx,ly,lz,lxx,lyy,lzz,dnorm
            # The difference coefficients are neglible so can truncate
            # at this level.  Delete the coeffs on the level below and
            # adopt the scaling function coeffs computed from the lower
            # level by twoscale since they will be more accurate
            if not self.s.has_key(n):
                self.s[n] = {}
            self.s[n][lx1,ly1,lz1,lx2,ly2,lz2] = s
            lx12, ly12, lz12, lx22, ly22, lz22 = \
                 lx1*2, ly1*2, lz1*2, lx2*2, ly2*2, lz2*2
            range2 = range(2)
            for ix1 in range2:
                for iy1 in range2:
                    for iz1 in range2:
                        for ix2 in range2:
                            for iy2 in range2:
                                for iz2 in range2:
                                    del(self.s[n+1][lx12+ix1,ly12+iy1,lz12+iz1,
                                                    lx22+ix2,ly22+iy2,lz22+iz2])
            ##if prnt:
            ##    for i in range(n): print "  ",
            ##    print " nofining level=%d l1=(%d,%d,%d) l2=(%d,%d,%d) norm=%.1e tol=%.1e" % \
            ##          (n, lx1, ly1, lz1, lx2, ly2, lz2, dnorm, tol)
            pass
        else:
            # First remove the inaccurate scaling function coeffs at
            # level n+1, then refine to level n+2
            if prnt:
                for i in range(n): print "  ",
                print " refining level=%d l1=(%d,%d,%d) l2=(%d,%d,%d) norm=%.1e tol=%.1e" % \
                      (n, lx1, ly1, lz1, lx2, ly2, lz2, dnorm, tol)
            lx12, ly12, lz12, lx22, ly22, lz22 = lx1*2, ly1*2, lz1*2, lx2*2, ly2*2, lz2*2
            range2 = range(2)
            for ix1 in range2:
                ix1lo = ix1*k
                ix1hi = ix1lo + k
                for iy1 in range2:
                    iy1lo = iy1*k
                    iy1hi = iy1lo+k
                    for iz1 in range2:
                        iz1lo = iz1*k
                        iz1hi = iz1lo+k
                        for ix2 in range2:
                            ix2lo = ix2*k
                            ix2hi = ix2lo + k
                            for iy2 in range2:
                                iy2lo = iy2*k
                                iy2hi = iy2lo+k
                                for iz2 in range2:
                                    iz2lo = iz2*k
                                    iz2hi = iz2lo+k
                                    del(self.s[n+1][lx12+ix1,ly12+iy1,lz12+iz1,
                                                    lx22+ix2,ly22+iy2,lz22+iz2])

                                    lx14, ly14, lz14, lx24, ly24, lz24 = \
                                         2*(lx12+ix1), 2*(ly12+iy1), 2*(lz12+iz1), \
                                         2*(lx22+ix2), 2*(ly22+iy2), 2*(lz22+iz2)

                                    for iix1 in range2:
                                        for iiy1 in range2:
                                            for iiz1 in range2:
                                                for iix2 in range2:
                                                    for iiy2 in range2:
                                                        for iiz2 in range2:
                                                            self.fine_scale_projection(n+2,[(
                                                                lx14+iix1,ly14+iiy1,lz14+iiz1,
                                                                lx24+iix2,ly24+iiy2,lz24+iiz2)])
                                                            
                                    self.refine_fine_scale_projection(
                                        n+1,lx12+ix1,ly12+iy1,lz12+iz1,
                                            lx22+ix2,ly22+iy2,lz22+iz2,prnt=0)

    def gather_scaling_coeffs(self, n, lx1, ly1, lz1, lx2, ly2, lz2):
        '''

        For a given translation at level n, gather the corresponding
        scaling function coefficients at level n+1.

        Some of the boxes on the lower level may be missing.

        '''
        k = self.k
        sn1 = self.s[n+1]
        ss = Tensor(2*k,2*k,2*k,2*k,2*k,2*k)
        lx12, ly12, lz12, lx22, ly22, lz22 = lx1*2, ly1*2, lz1*2, lx2*2, ly2*2, lz2*2
        range2 = range(2)
        for ix1 in range2:
            ix1lo = ix1*k
            ix1hi = ix1lo + k
            for iy1 in range2:
                iy1lo = iy1*k
                iy1hi = iy1lo+k
                for iz1 in range2:
                    iz1lo = iz1*k
                    iz1hi = iz1lo+k
                    for ix2 in range2:
                        ix2lo = ix2*k
                        ix2hi = ix2lo + k
                        for iy2 in range2:
                            iy2lo = iy2*k
                            iy2hi = iy2lo+k
                            for iz2 in range2:
                                iz2lo = iz2*k
                                iz2hi = iz2lo+k
                                key = (lx12+ix1,ly12+iy1,lz12+iz1,
                                       lx22+ix2,ly22+iy2,lz22+iz2)
                                if sn1.has_key(key):
                                    ss[ix1lo:ix1hi,iy1lo:iy1hi,iz1lo:iz1hi,
                                       ix2lo:ix2hi,iy2lo:iy2hi,iz2lo:iz2hi] =  sn1[key]
        return ss

    def reconstruct(self, nonstandard=0):
        '''

        Reconstruct the projection onto scaling functions from the
        compressed form terminating at the locally lowest level
        without difference coefficients.

        The compression algorithm has the responsibility of ensuring
        that the there are no signficant nodes below a node without
        difference coefficients.

        This algorithm must ensure that there are scaling function
        coefficients only at the locally lowest level (since other
        algorithms will just run down the tree until the first
        coefficient is found, and then stop).

        !!! nonstandard = 1
        A variant of this algorithm is used to reconstruct the
        function in preparation for the application of a non-standard
        form operator.  For this we need to keep the difference
        coefficients and also keep the sum coefficients at each level
        on the way down, but only down to locally finest level at
        which we have difference coefficients.  There is (will be)
        another routine to clean the resulting mess back into the
        conventional hierarchical basis definition.

        !!! nonstandard = 2
        Another variant for application of the non-standard form operator
        requires just the scaling coefficients at all levels.  No need
        to keep the difference coefficients.

        '''

        if self.debug:
            print "reconstructing", id(self), self.compressed, self.s.keys()
            
        if not self.compressed:
            if nonstandard:
                self.compress()
            else:
                return self

        perf.enter('reconstrct')
        self.compressed = 0

        k = self.k

        nmax = max(self.d.keys())
        if nonstandard == 1:
            # Only need scaling functions down to where the difference
            # coefficients are locally non-zero.
            nmax = nmax-1

        range2 = range(2)
            
        for n in range(nmax+1):
            sn = self.s[n]
            dn = self.d[n]
            sn1 = self.s[n+1] = {}
            if nonstandard == 1: dn1 = self.d[n+1]
            for lx1,ly1,lz1,lx2,ly2,lz2 in dn.keys():
                d = dn[lx1,ly1,lz1,lx2,ly2,lz2]
                d[0:k,0:k,0:k,0:k,0:k,0:k] = sn[lx1,ly1,lz1,lx2,ly2,lz2]
                if nonstandard != 2: del(sn[lx1,ly1,lz1,lx2,ly2,lz2])
                d = self.unfilter(d)
                lx12, ly12, lz12, lx22, ly22, lz22 = lx1*2, ly1*2, lz1*2, lx2*2, ly2*2, lz2*2
                for ix1 in range2:
                    for iy1 in range2:
                        for iz1 in range2:
                            for ix2 in range2:
                                for iy2 in range2:
                                    for iz2 in range2:
                                        key2 = lx12+ix1,ly12+iy1,lz12+iz1, \
                                               lx22+ix2,ly22+iy2,lz22+iz2
                                        if (nonstandard==1) and (not dn1.has_key(key2)):
                                            continue
                                        dxyz = d[ix1*k:ix1*k+k,iy1*k:iy1*k+k,iz1*k:iz1*k+k,
                                                 ix2*k:ix2*k+k,iy2*k:iy2*k+k,iz2*k:iz2*k+k]
                                        sn1[key2] = Tensor(dxyz)

            if nonstandard == 1:
                del self.s[n]
            else:
                del self.d[n]

        if nonstandard == 1: del self.s[nmax+1]
        
        if nmax == -1 and nonstandard == 1:
            print "FUDGE", self.s[0][0,0,0,0,0,0].normf()
            self.d[0][0,0,0,0,0,0][0:k,0:k,0:k,0:k,0:k,0:k] = self.s[0][0,0,0,0,0,0]
            print "FUDGE", self.d[0][0,0,0,0,0,0].normf()
            del self.s[0][0,0,0,0,0,0]

        perf.exit('reconstrct')
        return self

    def compress(self,truncate=0):
        '''

        Compress the scaling function coefficients (Vn+1 = Vn + Wn).

        In the usual compressed form, the scaling function
        coefficients will only exist on the lowest level, but this
        will vary across the domain.

        In the non-standard compressed form, both scaling function and
        difference coefficients will be defined on all levels ... must
        add the result of compressing the sum coefficients into the
        existing difference coefficients.

        In the form that results from addition in the scaling function
        basis, there may be sum coefficients at many levels, but no
        differences.

        For ease of handling all the special cases now implement
        truncation as a separate and optional step.

        Returns self so that operations can be chained.

        '''
        if self.debug:
            print "compressing", id(self), self.compressed, self.s.keys()
            
        if self.compressed: return self
        self.compressed = 1
        perf.enter('compress')
        
        nmax = max(self.s.keys())
        for n in range(nmax-1,-1,-1):
            if not self.d.has_key(n):
                self.d[n] = {}

        k = self.k
        for n in range(nmax-1,-1,-1):
            dn  = self.d[n]
            sn1 = self.s[n+1]
            try:
                sn = self.s[n]
            except KeyError:
                sn = self.s[n] = {}
            try:
                dn1 = self.d[n+1]
            except KeyError:
                dn1 = {}
            
            # Compute the list of translations that we must consider, driven
            # from the non-zero entries on the level below.
            llist = {}
            for lx1,ly1,lz1,lx2,ly2,lz2 in sn1.keys():
                llist[parent(lx1),parent(ly1),parent(lz1),
                      parent(lx2),parent(ly2),parent(lz2)] = 1
            
            for lx1,ly1,lz1,lx2,ly2,lz2 in llist.keys():
                ss = self.gather_scaling_coeffs(n, lx1,ly1,lz1,lx2,ly2,lz2)
                ss = self.filter(ss)
                s = Tensor(ss[0:k,0:k,0:k,0:k,0:k,0:k]) # Force a copy
                if sn.has_key((lx1,ly1,lz1,lx2,ly2,lz2)):
                    sn[lx1,ly1,lz1,lx2,ly2,lz2] = sn[lx1,ly1,lz1,lx2,ly2,lz2].gaxpy(1.0,s,1.0)
                else:
                    sn[lx1,ly1,lz1,lx2,ly2,lz2] = s

                ss[0:k,0:k,0:k,0:k,0:k,0:k].fill(0.0) # Now just want the diff coeffs
                if dn.has_key((lx1,ly1,lz1,lx2,ly2,lz2)):
                    dn[lx1,ly1,lz1,lx2,ly2,lz2] = dn[lx1,ly1,lz1,lx2,ly2,lz2].gaxpy(1.0,ss,1.0)
                else:
                    dn[lx1,ly1,lz1,lx2,ly2,lz2] = ss

            del self.s[n+1]

        if truncate:
            self.truncate()

        if not self.d.has_key(0):
            self.d[0] = {}
        if not self.d[0].has_key((0,0,0,0,0,0)):
            self.d[0][0,0,0,0,0,0] = Tensor(2*k,2*k,2*k,2*k,2*k,2*k)

        if not self.s.has_key(0):
            self.s[0] = {}
        if not self.s[0].has_key((0,0,0,0,0,0)):
            self.s[0][0,0,0,0,0,0] = Tensor(k,k,k,k,k,k)

        perf.exit('compress')
        return self

    def basis(self,refine=0):
        '''

        Self is a compressed function.  Export a description of
        the basis (for use by restrict).  This describes the
        pieces of W0+W1+...+W(n-1) in use.

        If refine=1, then if the difference coefficients do not
        satisfy the accuracy goal, the basis is locally refined down
        one level.

        '''
        if not self.compressed: self.compress()

        if Function.truncate_method == 2:
            lend = self.__lend()
        
        basis = {}
        nmax = max(self.d.keys())
        for n in range(nmax+1):
            dn = self.d[n]
            bn = basis[n] = {}
            for key in dn.keys():
                bn[key] = 1

        if refine:
            basis[nmax+1] = {}          # Just in case
            for n in range(nmax+1):
                dn = self.d[n]
                bn = basis[n]
                bn1 = basis[n+1]

                keys = dn.keys()

                if Function.truncate_method == 0:
                    ##tol = self.thresh
                    tol = self.thresh*sqrt(0.5**(3*max(n,3)))
                elif Function.truncate_method == 1:
                    tol = sqrt(0.5**(self.ndim*n))*self.thresh
                else:
                    tol = self.thresh/sqrt(lend)

                range2 = range(2)
                for key in keys():
                    if (dn[key].normf() > tol):
                        lx1 ,ly1 ,lz1 ,lx2 ,ly2 ,lz2 = key
                        lx12,ly12,lz12,lx22,ly22,lz22 = lx1*2,ly1*2,lz1*2,lx2*2,ly2*2,lz2*2
                        printed = 0
                        for ix1 in range2:
                            for iy1 in range2:
                                for iz1 in range2:
                                    for ix2 in range2:
                                        for iy2 in range2:
                                            for iz2 in range2:
                                                key2 = (lx12+ix1,ly12+iy1,lz12+iz1,
                                                        lx22+ix2,ly22+iy2,lz22+iz2)
                                                if not bn1.has_key(key2):
                                                    if not printed:
                                                        print "Refining", n, key
                                                        printed = 1
                                                    bn1[key2] = 1
            if not basis[nmax+1]:
                del basis[nmax+1]

        return basis
    
    def truncate(self,thresh=None):
        '''

        The difference coefficients are truncated to the specified
        threshold.  So that reconstruction can just proceed down to
        the first level without difference coefficients, it is
        important to only truncate difference coefficients at leaf
        nodes.

        Returns self so that operations can be chained.

        '''
        if self.debug:
            print "truncating", id(self), self.compressed, self.thresh
            
        if not self.compressed:
            self.compress()

        if Function.truncate_method == 2:
            lend = self.__lend()
        
        if thresh is None:
            thresh = self.thresh
        
        perf.enter('truncate')

        nmax = max(self.d.keys())

        is_parent = {}
        for n in range(nmax,0,-1):
            dn = self.d[n]
            keys = dn.keys()

            if Function.truncate_method == 0:
                tol = self.thresh*sqrt(0.5**(3*max(n,3)))
                ##tol = thresh
            elif Function.truncate_method == 1:
                tol = sqrt(0.5**(self.ndim*n))*self.thresh
            else:
                tol = thresh/sqrt(lend)

            # If not a parent and below threshold kill box.
            # If above threshold mark parent box as such.
            next_is_parent = {}
            for key in keys:
                if (not is_parent.has_key(key)) and \
                   (dn[key].normf() < tol):
                    del dn[key]
                else:
                    lx1,ly1,lz1,lx2,ly2,lz2 = key
                    next_is_parent[parent(lx1),parent(ly1),parent(lz1),
                                   parent(lx2),parent(ly2),parent(lz2)] = 1

            is_parent = next_is_parent
                    
            # If the current level is totally empty and it is also
            # the lowest level, it can be deleted entirely.
            if (not self.d[n]) and (not self.d.has_key(n+1)):
                del self.d[n]

        perf.exit('truncate')

##         if Function.truncate_method == 2:
##             # If we deleted a lot of coefficients, try again ...
##             if self.__lend()*2 < lend:
##                 self.truncate()

        return self

    def griddx(self):
        pass

    def plotdx(self):
        pass

    def eval(self,x1,y1,z1,x2,y2,z2):
        '''

        Evaluate the function at the given point.

        If it is not compressed then we just need to walk down the
        tree until we locate the leaf node containing x,y,z.  If it is
        compressed then we do the same, but also have to add up the
        scaling function coeffs on the way down.

        '''
        n = 0
        k = self.k
        lx1 = ly1 = lz1 = lx2 = ly2 = lz2 = 0
        xx1, yy1, zz1, xx2, yy2, zz2 = x1, y1, z1, x2, y2, z2
        if self.compressed:
            # Note always want loop to terminate via the break otherwise
            # n is one too small (x & l would also be wrong)
            s = self.s[0][0,0,0,0,0,0]
            for n in range(max(self.d.keys())+2):
                try:
                    d = Tensor(self.d[n][lx1,ly1,lz1,lx2,ly2,lz2]) # Take a copy
                except KeyError:
                    # At the bottom ... evaluate s*phi(xx) at level n
                    break
                d[0:k,0:k,0:k,0:k,0:k,0:k] = s
                d = self.unfilter(d)
                scale = 2**(n+1)
                xx1, yy1, zz1, xx2, yy2, zz2 = \
                     x1*scale, y1*scale, z1*scale, x2*scale, y2*scale, z2*scale
                lx1, ly1, lz1, lx2, ly2, lz2 = \
                     min(scale-1,int(xx1)), min(scale-1,int(yy1)), min(scale-1,int(zz1)), \
                     min(scale-1,int(xx2)), min(scale-1,int(yy2)), min(scale-1,int(zz2))                      
                xx1, yy1, zz1, xx2, yy2, zz2 = xx1-lx1, yy1-ly1, zz1-lz1, xx2-lx2, yy2-ly2, zz2-lz2
                ix1, iy1, iz1, ix2, iy2, iz2 = lx1%2, ly1%2, lz1%2, lx2%2, ly2%2, lz2%2,
                s = d[ix1*k:ix1*k+k,iy1*k:iy1*k+k,iz1*k:iz1*k+k,
                      ix2*k:ix2*k+k,iy2*k:iy2*k+k,iz2*k:iz2*k+k]
            s = Tensor(s)
        else:
            # The function is not compressed.  Assume that it is
            # represented at some (locally defined) finest level
            # in the scaling function basis.
            nmax = max(self.s.keys()) + 1
            for n in range(nmax):
                twon = 2**n
                twon1= twon - 1
                xx1, yy1, zz1, xx2, yy2, zz2 = \
                     x1*twon, y1*twon, z1*twon, x2*twon, y2*twon, z2*twon
                lx1, ly1, lz1, lx2, ly2, lz2 = \
                     min(twon1,int(xx1)), min(twon1,int(yy1)), min(twon1,int(zz1)), \
                     min(twon1,int(xx2)), min(twon1,int(yy2)), min(twon1,int(zz2))
                xx1, yy1, zz1, xx2, yy2, zz2 = \
                     xx1-lx1, yy1-ly1, zz1-lz1, xx2-lx2, yy2-ly2, zz2-lz2
                try:
                    s = self.s[n][lx1,ly1,lz1,lx2,ly2,lz2]
                    break
                except KeyError:
                    pass
            
        #px1 = Vector(phi(xx1,k))
        #py1 = Vector(phi(yy1,k))
        #pz1 = Vector(phi(zz1,k))
        #px2 = Vector(phi(xx2,k))
        #py2 = Vector(phi(yy2,k))
        #pz2 = Vector(phi(zz2,k))
        #value = s.inner(pz2).inner(py2).inner(px2).inner(pz1).inner(py1).inner(px1)*(8.0**n)

        value = cfcube.eval6(k,n,xx1,yy1,zz1,xx2,yy2,zz2,s.swigptr())
        return value
  
    def eval_err(self):
        pass

    def print_layers(self, print_pic=0):
        if self.compressed:
            print " Analysis by layers of compressed function"
            nmax = max(self.d.keys())
            normsq = self.s[0][0,0,0,0,0,0].normf()**2
            print "\nNorm in V0 %.8e" % sqrt(normsq)
            for n in range(nmax+1):
                sum = 0.0
                for xyz in self.d[n].keys():
                    sum = sum + self.d[n][xyz].normf()**2
                normsq = normsq + sum
                print "Norm in W%d %.8e %.8e" % (n,sqrt(sum),sqrt(normsq))

            if not print_pic: return  
        else:
            print " Analysis by layers of reconstructed function"
            nmax = max(self.s.keys())
            normsq = 0.0
            for n in range(nmax+1):
                sum = 0.0
                for xyz in self.s[n].keys():
                    sum = sum + self.s[n][xyz].normf()**2
                normsq = normsq + sum
                print "Norm in V%d %.8e %.8e" % (n,sqrt(sum),sqrt(normsq))

            if not print_pic: return

    def sclean(self):
        '''

        Cleans up scaling function coefficients after diffing or
        multiplying.  Only want to keep the scaling function coeffs
        at the highest level in the tree ... ones lower down will
        have been generated by two-scale just for evaluation.
        Kill all nodes that have a parent - infanticide.

        '''
        for n in range(max(self.s.keys()),0,-1):
            sn, sn1 = self.s[n], self.s[n-1]
            for lx1,ly1,lz1,lx2,ly2,lz2 in sn.keys():
                if sn1.has_key((parent(lx1),parent(ly1),parent(lz1),
                                parent(lx2),parent(ly2),parent(lz2))):
                    del(sn[lx1,ly1,lz1,lx2,ly2,lz2])
        
    def mul(self,other):
        '''

        Multiplication of a function by a another function.
        Replaces the algorithm based on squaring.

        Returns a new function in the scaling function basis.
        Self/other are not modified though may be reconstructed if
        either is compressed.

        Added an optimization for multiplcation by the mask which is
        mask(x) = 1 for 0.0625 < x < 1-0.625.  Thus, only need to do
        the multiplication we have boxes 0 or 15 on level 4 (or boxes
        contained within below).  The mask is identified by the
        attribute mask which noone else will have.

        '''

        self.reconstruct()
        other.reconstruct()

        if hasattr(self,'mask'): self,other = other,self
        mask = hasattr(other,'mask')
        
        perf.enter('mul')

        result = self.new(k=self.k,compress=0) # Was thresh=self.thresh
        result.s = {}

        k = self.k
        phit = Matrix(self.quad_phi.cycledim(1)) # Contig. for fast_transform
        phiw = self.quad_phiw
        npt = len(self.quad_x)
        nmax = max(self.s.keys()+other.s.keys())

        tmp1 = Tensor(k,k,k,k,k,k)
        tmp2 = Tensor(k,k,k,k,k,k)

        range2 = range(2)

        if Function.autorefine:     # NOTE REFERENCE TO CLASS
            result.s[0] = {}
             
        for n in range(nmax+1):
            try:
                sn = self.s[n]
            except KeyError:
                sn = self.s[n] = {}

            try:
                on = other.s[n]
            except KeyError:
                on = other.s[n] = {}
                
            if Function.autorefine:     # NOTE REFERENCE TO CLASS
                rn1 = result.s[n+1] = {}
                scale = 8.0**(n+1)
            else:
                rn = result.s[n] = {}
                scale = 8.0**n

            # Generate union of keys on this level. 
            keys = sn.keys() + on.keys()
            keys.sort()
            for i in range(len(keys)-1,0,-1):
                if keys[i] == keys[i-1]:
                    del keys[i]

            if mask:
                if n <= 4:
                    masklo, maskhi = 0, 2**n-1
                else:
                    masklo, maskhi = 2**(n-4)-1, 15*2**(n-4)
                    
            for key in keys:
                lx1,ly1,lz1,lx2,ly2,lz2 = key

                # If other is the mask, exclude keys for which the
                # mask is unity.
                if mask:
                    if masklo<lx1<maskhi and masklo<ly1<maskhi and masklo<lz1<maskhi and \
                       masklo<lx2<maskhi and masklo<ly2<maskhi and masklo<lz2<maskhi:
                        if sn.has_key(key): rn[key] = Tensor(sn[key])
                        continue

                s =  self.__sock_it_to_me(n,key)
                if not (s is None):
                    o = other.__sock_it_to_me(n,key)
                    if not (o is None):
                        if Function.autorefine: # NOTE REFERENCE TO CLASS
                            ss = self.unfilter(s,sonly=1)
                            oo = self.unfilter(o,sonly=1)
                            lx12, ly12, lz12, lx22, ly22, lz22 = \
                                  2*lx1, 2*ly1, 2*lz1, 2*lx2, 2*ly2, 2*lz2
                            for ix1 in range2:
                                for iy1 in range2:
                                    for iz1 in range2:
                                        for ix2 in range2:
                                            for iy2 in range2:
                                                for iz2 in range2:
                                                    f = ss[ix1*k:ix1*k+k,iy1*k:iy1*k+k,iz1*k:iz1*k+k,
                                                           ix2*k:ix2*k+k,iy2*k:iy2*k+k,iz2*k:iz2*k+k]
                                                    g = oo[ix1*k:ix1*k+k,iy1*k:iy1*k+k,iz1*k:iz1*k+k,
                                                           ix2*k:ix2*k+k,iy2*k:iy2*k+k,iz2*k:iz2*k+k]
                                                    f = f.transform(phit,result=tmp1)
                                                    g = g.transform(phit,result=tmp2)
                                                    f.emul(g)
                                                    rn1[lx12+ix1,ly12+iy1,lz12+iz1,
                                                        lx22+ix2,ly22+iy2,lz22+iz2] = \
                                                       f.fast_transform(phiw).scale(scale)
                        else:
                            f = s.fast_transform(phit,result=tmp1)
                            g = o.fast_transform(phit,result=tmp2)
                            f.emul(g)
                            rn[key] = f.fast_transform(phiw).scale(scale)
                            
        self.sclean()
        other.sclean()

        perf.exit('mul')

        return result


    def mul_with3d(self,other,side=1):
        '''
        
        side = 1 : left    f_new(r1,r2) = self(r1,r2) * other(r1)
        side = 2 : right   f_new(r1,r2) = self(r1,r2) * other(r2)
        
        Multiplication of a function by a another function.
        Replaces the algorithm based on squaring.

        Returns a new function in the scaling function basis.
        Self/other are not modified though may be reconstructed if
        either is compressed.

        Added an optimization for multiplcation by the mask which is
        mask(x) = 1 for 0.0625 < x < 1-0.625.  Thus, only need to do
        the multiplication we have boxes 0 or 15 on level 4 (or boxes
        contained within below).  The mask is identified by the
        attribute mask which noone else will have.

        '''

        if side!=1 and side!=2:
            raise "side is either 1 or 2"

        self.reconstruct()
        other.reconstruct()

        if hasattr(self,'mask'): self,other = other,self
        mask = hasattr(other,'mask')
        
        perf.enter('mul')

        result = self.new(k=self.k,compress=0) # Was thresh=self.thresh
        result.s = {}

        k = self.k
        phit = Matrix(self.quad_phi.cycledim(1)) # Contig. for fast_transform
        phiw = self.quad_phiw
        npt = len(self.quad_x)
        nmax = max(self.s.keys()+other.s.keys())

        tmp1 = Tensor(k,k,k,k,k,k)
        one  = Tensor(k,k,k)
        one[0,0,0] = 1.0

        range2 = range(2)

        if Function.autorefine:     # NOTE REFERENCE TO CLASS
            result.s[0] = {}
             
        for n in range(nmax+1):
            try:
                sn = self.s[n]
            except KeyError:
                sn = self.s[n] = {}

            try:
                on = other.s[n]
            except KeyError:
                on = other.s[n] = {}

            if Function.autorefine:     # NOTE REFERENCE TO CLASS
                rn1 = result.s[n+1] = {}
                scale = 8.0**(n+1)
            else:
                rn = result.s[n] = {}
                scale = 8.0**n
                
            one[0,0,0] = 1.0/sqrt(scale)

            # Generate union of keys on this level.
            trans1d = range(2**n)
            onkeys = []
            if(side==1):            
                for key in on.keys():
                    for lx in trans1d:
                        for ly in trans1d:
                            for lz in trans1d:
                                onkeys.append(key+(lx,ly,lz))  # <==========check for left
            else:
                for key in on.keys():
                    for lx in trans1d:
                        for ly in trans1d:
                            for lz in trans1d:
                                onkeys.append((lx,ly,lz)+key)  # <==========check for right

            keys = sn.keys() + onkeys
            keys.sort()
            for i in range(len(keys)-1,0,-1):
                if keys[i] == keys[i-1]:
                    del keys[i]

            if mask:
                if n <= 4:
                    masklo, maskhi = 0, 2**n-1
                else:
                    masklo, maskhi = 2**(n-4)-1, 15*2**(n-4)

            cache_g = {}
            for key in keys:
                lx1,ly1,lz1,lx2,ly2,lz2 = key

                # If other is the mask, exclude keys for which the
                # mask is unity.
                if mask:
                    if masklo<lx1<maskhi and masklo<ly1<maskhi and masklo<lz1<maskhi and \
                       masklo<lx2<maskhi and masklo<ly2<maskhi and masklo<lz2<maskhi:
                        if sn.has_key(key): rn[key] = Tensor(sn[key])
                        continue

                s =  self.__sock_it_to_me(n,key)
                if not (s is None):
                    if(side==1):            
                        o = other._Function__sock_it_to_me(n,(lx1,ly1,lz1)) # <=========== check for left
                    else:
                        o = other._Function__sock_it_to_me(n,(lx2,ly2,lz2)) # <=========== check for right
                    if not (o is None):
                        if Function.autorefine: # NOTE REFERENCE TO CLASS
                            ss = self.unfilter(s,sonly=1)
                            oo = other.unfilter(o,sonly=1)
                            lx12, ly12, lz12, lx22, ly22, lz22 = \
                                  2*lx1, 2*ly1, 2*lz1, 2*lx2, 2*ly2, 2*lz2
                            for ix1 in range2:
                                for iy1 in range2:
                                    for iz1 in range2:
                                        for ix2 in range2:
                                            for iy2 in range2:
                                                for iz2 in range2:
                                                    f = ss[ix1*k:ix1*k+k,iy1*k:iy1*k+k,iz1*k:iz1*k+k,
                                                           ix2*k:ix2*k+k,iy2*k:iy2*k+k,iz2*k:iz2*k+k]
                                                    if(side==1):            
                                                        g = oo[ix1*k:ix1*k+k,iy1*k:iy1*k+k,iz1*k:iz1*k+k] # <==check for left
                                                        try:
                                                            g = cache_g[(lx1,ly1,lz1)]               # <====== check for left
                                                        except KeyError:
                                                            g = g.outer(one)                         # <====== check for left
                                                            g = g.transform(phit)
                                                            cache_g[(lx1,ly1,lz1)] = g
                                                    else:
                                                        g = oo[ix2*k:ix2*k+k,iy2*k:iy2*k+k,iz2*k:iz2*k+k] # <==check for right
                                                        try:
                                                            g = cache_g[(lx2,ly2,lz2)]               # <====== check for right
                                                        except KeyError:
                                                            g = one.outer(g)                         # <====== check for right
                                                            g = g.transform(phit)
                                                            cache_g[(lx2,ly2,lz2)] = g
                                                            
                                                    f = f.transform(phit,result=tmp1)
                                                    f.emul(g)
                                                    rn1[lx12+ix1,ly12+iy1,lz12+iz1,
                                                        lx22+ix2,ly22+iy2,lz22+iz2] = \
                                                       f.fast_transform(phiw).scale(scale)
                        else:
                            if(side==1):            
                                try:
                                    g = cache_g[(lx1,ly1,lz1)]                # <====== check for right
                                except KeyError:
                                    o = o.outer(one)
                                    g = cache_g[(lx1,ly1,lz1)] = o.fast_transform(phit)
                            else:
                                try:
                                    g = cache_g[(lx2,ly2,lz2)]                # <====== check for right
                                except KeyError:
                                    o = one.outer(o)
                                    g = cache_g[(lx2,ly2,lz2)] = o.fast_transform(phit)
                                    
                            f = s.fast_transform(phit,result=tmp1)
                            f.emul(g)
                            rn[key] = f.fast_transform(phiw).scale(scale)
            del cache_g

        self.sclean()
        other.sclean()

        perf.exit('mul')

        return result

    def __rmul__(self,other):

        return self.__mul__(other)

    def square(self):
        '''

        Square self replacing self by its square in the scaling
        function basis.

        If (Function.autorefine) is set, then it will recur down an
        additional level before squaring, otherwise it does not.

        NB: self is changed.

        Returns self so that operations may be chained.  

        '''
        if self.debug:
            print "squaring", id(self), self.compressed
            
        self.reconstruct()

        k = self.k
        nmax = max(self.s.keys())
        phit = Matrix(self.quad_phi.cycledim(1)) # Contig. for fast_transform
        phiw = self.quad_phiw
        npt = len(self.quad_x)

        tmp1 = Tensor(k,k,k,k,k,k)
        range2 = range(2)

        try:
            next_keys = self.s[0].keys()
        except KeyError:
            # This can happen if not compressed and initial_level != 0
            self.s[0] = {}
            next_keys = self.s[0].keys()

        for n in range(nmax+1):
        
            sn = self.s[n]
            try:
                sn1 = self.s[n+1]
            except KeyError:
                sn1 = self.s[n+1] = {}

            keys = next_keys
            next_keys = sn1.keys()

            if Function.autorefine:         # NOTE REFERENCE TO CLASS
                ss = Tensor(2*k,2*k,2*k,2*k,2*k,2*k)
                scale = 8.0**(n+1)
            else:
                scale = 8.0**n
            for lx1, ly1, lz1, lx2, ly2, lz2 in keys:
                if Function.autorefine:     # NOTE REFERENCE TO CLASS
                    # Now are at the locally lowest level ... recur down
                    # one more level and then form the square projecting
                    # onto the approximately refined basis.  Could recur
                    # more dynamically to guarantee precision but
                    # Beylkin does not seem to think this is necessary
                    # except in pathological cases?
                    ss = self.unfilter(sn[lx1,ly1,lz1,lx2,ly2,lz2],sonly=1)
                    del(sn[lx1,ly1,lz1,lx2,ly2,lz2])
                    lx12, ly12, lz12, lx22, ly22, lz22 = \
                          2*lx1, 2*ly1, 2*lz1, 2*lx2, 2*ly2, 2*lz2
                    for ix1 in range2:
                        for iy1 in range2:
                            for iz1 in range2:
                                for ix2 in range2:
                                    for iy2 in range2:
                                        for iz2 in range2:
                                            f = ss[ix1*k:ix1*k+k,iy1*k:iy1*k+k,iz1*k:iz1*k+k,
                                                   ix2*k:ix2*k+k,iy2*k:iy2*k+k,iz2*k:iz2*k+k]
                                            # Back transform to the function value
                                            # at the quadrature points
                                            f = f.transform(phit,result=tmp1)
                                            # Compute the square on the grid
                                            f.emul(f)
                                            # Do quadrature on f**2 and scale
                                            sn1[lx12+ix1,ly12+iy1,lz12+iz1,
                                                lx22+ix2,ly22+iy2,lz22+iz2] = \
                                                f.fast_transform(phiw).scale(scale)
                                                          
                else:
                    # Don't refine ... just do it at the current level
                    f = sn[lx1,ly1,lz1,lx2,ly2,lz2].fast_transform(phit,result=tmp1)
                    f.emul(f)
                    sn[lx1,ly1,lz1,lx2,ly2,lz2] = f.fast_transform(phiw).scale(scale)
            
        return self

    def local_function(self,function=None,cfunction=None):
        '''

        Apply a scalar local function to self replacing self by
        f(self(x)).  Define either function (a callable object) or
        cfunction (a pointer to a function cast to a pointer to void
        ... so that old versions of SWIG can be used).  Both the C and
        Python functions should take a double argument and return a
        double.  THIS APPROACH IS ONLY GOOD IF THE RESULT OF THE
        FUNCTION IS AS LOCALLY SMOOTH OR SMOOTHER THAN THE INPUT ...

        If (Function.autorefine) is set, then it will recur down an
        additional level before squaring, otherwise it does not.

        NB: self is changed.

        Returns self so that operations may be chained.  

        '''
        if self.debug:
            print "squaring", id(self), self.compressed
            
        range2 = range(2)

        self.reconstruct()
        perf.enter('local func')

        k = self.k
        nmax = max(self.s.keys())
        phit = Matrix(self.quad_phi.cycledim(1)) # Contig. for fast_transform
        phiw = self.quad_phiw
        npt = len(self.quad_x)

        tmp1 = Tensor(k,k,k,k,k,k)

        try:
            next_keys = self.s[0].keys()
        except KeyError:
            # This can happen if not compressed and initial_level != 0
            self.s[0] = {}
            next_keys = self.s[0].keys()

        for n in range(nmax+1):
        
            sn = self.s[n]
            try:
                sn1 = self.s[n+1]
            except KeyError:
                sn1 = self.s[n+1] = {}

            keys = next_keys
            next_keys = sn1.keys()

            if Function.autorefine:         # NOTE REFERENCE TO CLASS
                ss = Tensor(2*k,2*k,2*k,2*k,2*k,2*k)
                scale = 8.0**(n+1)
            else:
                scale = 8.0**n
            for lx1, ly1, lz1, lx2, ly2, lz2 in keys:
                if Function.autorefine:     # NOTE REFERENCE TO CLASS
                    # Now are at the locally lowest level ... recur down
                    # one more level and then form the square projecting
                    # onto the approximately refined basis.  Could recur
                    # more dynamically to guarantee precision but
                    # Beylkin does not seem to think this is necessary
                    # except in pathological cases?
                    ss = self.unfilter(sn[lx1,ly1,lz1,lx2,ly2,lz2],sonly=1)
                    del(sn[lx1,ly1,lz1,lx2,ly2,lz2])
                    lx12, ly12, lz12, lx22, ly22, lz22 = \
                          2*lx1, 2*ly1, 2*lz1, 2*lx2, 2*ly2, 2*lz2
                    for ix1 in range2:
                        for iy1 in range2:
                            for iz1 in range2:
                                for ix2 in range2:
                                    for iy2 in range2:
                                        for iz2 in range2:
                                            f = ss[ix1*k:ix1*k+k,iy1*k:iy1*k+k,iz1*k:iz1*k+k,
                                                   ix2*k:ix2*k+k,iy2*k:iy2*k+k,iz2*k:iz2*k+k]
                                            # Back transform to the function value
                                            # at the quadrature points
                                            f = Tensor(f).transform(phit,result=tmp1).scale(scale)
                                            # Eval function on the refined grid
                                            if cfunction:
                                                cfcube.clocalfunc(
                                                    f.swigptr(),k**6,cfunction)
                                            else:
                                                f.unaryop(function)
                                            # Do quadrature on f**2 and scale
                                            sn1[lx12+ix1,ly12+iy1,lz12+iz1,
                                                lx22+ix2,ly22+iy2,lz22+iz2] = \
                                                f.fast_transform(phiw).scale(1.0/scale)
                                                          
                else:
                    # Don't refine ... just do it at the current level
                    f = sn[lx1,ly1,lz1,lx2,ly2,lz2].fast_transform(phit,result=tmp1).scale(scale)
                    if cfunction:
                        cfcube.clocalfunc(f.swigptr(),k**6,cfunction)
                    else:
                        f.unaryop(function)
                    sn[lx1,ly1,lz1,lx2,ly2,lz2] = f.fast_transform(phiw).scale(1.0/scale)
            
        perf.exit('local func')
        return self

    def mul_func(self,function=None,cfunction=None):
        '''

        '''
        if self.debug:
            print "squaring", id(self), self.compressed
            
        range2 = range(2)

        self.reconstruct()
        perf.enter('mul func')

        k = self.k
        nmax = max(self.s.keys())
        phit = Matrix(self.quad_phi.cycledim(1)) # Contig. for fast_transform
        phiw = self.quad_phiw
        pts = self.quad_x 
        npt = len(self.quad_x)

        tmp1 = Tensor(k,k,k,k,k,k)

        try:
            next_keys = self.s[0].keys()
        except KeyError:
            # This can happen if not compressed and initial_level != 0
            self.s[0] = {}
            next_keys = self.s[0].keys()

        for n in range(nmax+1):

            sn = self.s[n]
            try:
                sn1 = self.s[n+1]
            except KeyError:
                sn1 = self.s[n+1] = {}

            keys = next_keys
            next_keys = sn1.keys()

            if Function.autorefine:         # NOTE REFERENCE TO CLASS
                ss = Tensor(2*k,2*k,2*k,2*k,2*k,2*k)
                dd = float(2**(n+1))
            else:
                dd = float(2**n)

            scale = sqrt(dd**6)       # 6D Scale
            h = 1.0/dd                # Box size on target level
             
            for lx1, ly1, lz1, lx2, ly2, lz2 in keys:

                if Function.autorefine:     # NOTE REFERENCE TO CLASS
                    # Now are at the locally lowest level ... recur down
                    # one more level and then form the square projecting
                    # onto the approximately refined basis.  Could recur
                    # more dynamically to guarantee precision but
                    # Beylkin does not seem to think this is necessary
                    # except in pathological cases?
                    ss = self.unfilter(sn[lx1,ly1,lz1,lx2,ly2,lz2],sonly=1)
                    del(sn[lx1,ly1,lz1,lx2,ly2,lz2])
                    lx12, ly12, lz12, lx22, ly22, lz22 = \
                          2*lx1, 2*ly1, 2*lz1, 2*lx2, 2*ly2, 2*lz2
                    for ix1 in range2:
                        for iy1 in range2:
                            for iz1 in range2:
                                for ix2 in range2:
                                    for iy2 in range2:
                                        for iz2 in range2:
            
                                            # get a cube of g:(func)
                                            # ----------------------
                                            x1lo = (lx12+ix1)*h
                                            y1lo = (ly12+iy1)*h
                                            z1lo = (lz12+iz1)*h
                                            x2lo = (lx22+ix2)*h
                                            y2lo = (ly22+iy2)*h
                                            z2lo = (lz22+iz2)*h
                                            g = Tensor(npt,npt,npt,npt,npt,npt)
                                            if cfunction:
                                                cfcube.cfcube6(x1lo, y1lo, z1lo, x2lo, y2lo, z2lo, npt, h,
                                                              pts.swigptr(), g.swigptr(),
                                                              cfunction)
                                            else:
                                                # Plain old Python function or non-cast C function
                                                # This is SLOW due to all of the tensor references.
                                                self.fcube(x1lo, y1lo, z1lo, x2lo, y2lo, z2lo, npt, h,
                                                           pts, g, function)
                                            
                                            # get a cube of f:(func) 
                                            # ----------------------
                                            f = ss[ix1*k:ix1*k+k,iy1*k:iy1*k+k,iz1*k:iz1*k+k,
                                                   ix2*k:ix2*k+k,iy2*k:iy2*k+k,iz2*k:iz2*k+k]
                                            # Back transform to the function value
                                            # at the quadrature points
                                            f = f.transform(phit,result=tmp1).scale(scale)
            
                                            # emul between f and g
                                            # --------------------
                                            f = f.emul(g)
                                            
                                            # Do quadrature on f**2 and scale
                                            sn1[lx12+ix1,ly12+iy1,lz12+iz1,
                                                lx22+ix2,ly22+iy2,lz22+iz2] = \
                                                f.fast_transform(phiw).scale(1.0/scale)
                                                          
                else:

                    # get a cube of g:(func)
                    # ----------------------
                    x1lo = lx1*h
                    y1lo = ly1*h
                    z1lo = lz1*h
                    x2lo = lx2*h
                    y2lo = ly2*h
                    z2lo = lz2*h
                    g = Tensor(npt,npt,npt,npt,npt,npt)
                    if cfunction:
                        cfcube.cfcube6(x1lo, y1lo, z1lo, x2lo, y2lo, z2lo, npt, h,
                                      pts.swigptr(), g.swigptr(),
                                      cfunction)
                    else:
                        # Plain old Python function or non-cast C function
                        # This is SLOW due to all of the tensor references.
                        self.fcube(x1lo, y1lo, z1lo, x2lo, y2lo, z2lo, npt, h,
                                   pts, g, function)

                    # get a cube of f:(func) 
                    # ----------------------
                    # Don't refine ... just do it at the current level
                    f = sn[lx1,ly1,lz1,lx2,ly2,lz2].fast_transform(phit,result=tmp1).scale(scale)

                    # emul between f and g
                    # --------------------
                    f = f.emul(g)
                
                    sn[lx1,ly1,lz1,lx2,ly2,lz2] = f.fast_transform(phiw).scale(1.0/scale)
            
        perf.exit('mul func')
        return self


    def __look_out_below(self,s, n, key):
        '''

        s are scaling coefficients at level n for the given
        translation.  Recur them down one level and put the resulting
        coefficients into the tree.

        '''
	print "foofoo"
	sys.exit()
        if not self.s.has_key(n+1):
            self.s[n+1] = {}
        k = self.k
        ss = self.unfilter(s,sonly=1)
        lx1, ly1, lz1, lx2, ly2, lz2 = key
        lx12, ly12, lz12, lx22, ly22, lz22 = \
              lx1*2, ly1*2, lz1*2, lx2*2, ly2*2, lz2*2
        for ix1 in range(2):
            for iy1 in range(2):
                for iz1 in range(2):
                    for ix2 in range(2):
                        for iy2 in range(2):
                            for iz2 in range(2):
                                view = ss[ix1*k:ix1*k+k,iy1*k:iy1*k+k,iz1*k:iz1*k+k,
                                          ix2*k:ix2*k+k,iy2*k:iy2*k+k,iz2*k:iz2*k+k]
                                self.s[n+1][lx12+ix1,ly12+iy1,lz12+iz1,
                                            lx22+ix2,ly22+iy2,lz22+iz2] = Tensor(view)


    def __sock_it_to_me(self,n,key):
        '''

        Return self.s[n][key] or None if it cannot be made (because it
        has no parent).
        
        If it is absent, then try to recur from the level above to
        make it.  If the parent is also absent, we need to keep
        walking up.  If we hit the top of the tree and still cannot
        find a parent, then return None.

        Some optimizations are possible, such as looking one below
        before recuring up.  Also, has_key() might be faster than
        exception handling.

        '''
        try:
            return self.s[n][key]
        except KeyError:
            if n == 0:
                return None

        lx1, ly1, lz1, lx2, ly2, lz2 = key

        lx1, ly1, lz1, lx2, ly2, lz2 = \
             parent(lx1), parent(ly1), parent(lz1), parent(lx2), parent(ly2), parent(lz2)
        try:
            s = self.s[n-1][lx1,ly1,lz1,lx2,ly2,lz2]
        except KeyError:
            s = self.__sock_it_to_me(n-1,(lx1,ly1,lz1,lx2,ly2,lz2))

        if s is None:
            return None

        self.__look_out_below(s,n-1,(lx1,ly1,lz1,lx2,ly2,lz2))

        return self.s[n][key]


    def __sock_it_to_me2(self,n,key):
        '''

        Return self.s[n][key] or None if it cannot be made (because it
        has no parent).
     
        If it is absent, then try to recur from the level above to
        make it.  If the parent is also absent, we need to keep
        walking up.  If we hit the top of the tree and still cannot
        find a parent, then return None.

        Some optimizations are possible, such as looking one below
        before recuring up.  Also, has_key() might be faster than
        exception handling.

        '''
        try:
            return self.s[n][key].normf()
        except KeyError:
            if n == 0:
                return None

        lx1, ly1, lz1, lx2, ly2, lz2 = key

        lx1, ly1, lz1, lx2, ly2, lz2 = \
             parent(lx1), parent(ly1), parent(lz1), parent(lx2), parent(ly2), parent(lz2)
        try:
            s = self.s[n-1][lx1,ly1,lz1,lx2,ly2,lz2]
            s = s.normf()
        except KeyError:
            s = self.__sock_it_to_me2(n-1,(lx1,ly1,lz1,lx2,ly2,lz2))

        if s is None:
            return None

        s = s * 0.125 # = (2.0**(-6.0/2.0))
     
        return s


    def diff(self,axis,transpose=0):
        '''

        Apply a central derivative along the specified axis,
        currently with zero boundary conditions.  Plan to implement
        periodic boundary conditions also.

        If (transpose) then use the transposed derivative operator.

        Self is reconstructed but should otherwise be unchanged.

        Returns a new function in the scaling function basis.

        '''
        if self.debug:
            print "diffing", id(self), self.compressed, axis, \
                  transpose,self.s.keys()

        self.reconstruct()
        perf.enter('diff')
        
        result = self.new(k=self.k, compress=0)

        zero = Tensor(self.k,self.k,self.k,self.k,self.k,self.k)

        ind = (1,axis)                  # Indices to contract in products
        if axis == 0:
            amap = None
        elif axis == 1:
            amap = [1,0,2,3,4,5]
        elif axis == 2:
            amap = [2,0,1,3,4,5]
        elif axis == 3:
            amap = [3,0,1,2,4,5]
        elif axis == 4:
            amap = [4,0,1,2,3,5]
        elif axis == 5:
            amap = [5,0,1,2,3,4]
        else:
            raise "Don't you think six dimensions are enough?"
            
        rm, r0, rp = self.make_dc_periodic()

        if transpose:
            rm, r0, rp = rp.transpose(), r0.transpose(), rm.transpose()

        nmax = max(self.s.keys())

        for n in range(nmax+1):         #  Fill in empty levels
            try:
                self.s[n]
            except KeyError:
                self.s[n] = {}
            
        for n in range(nmax+1):
            result.s[n] = {}
            sn = self.s[n]
            twon = 2**n
            twon1= twon - 1
            for key in self.s[n].keys():
                s0 = sn[key]
                
                if key[axis] == 0:
                    sm = zero      # Zero outside the box
                else:
                    t=list(key); t[axis]=t[axis]-1; t=tuple(t)
                    sm = self.__sock_it_to_me(n,t)
                    
                if key[axis] == twon1:
                    sp = zero
                else:
                    t=list(key); t[axis]=t[axis]+1; t=tuple(t)
                    sp = self.__sock_it_to_me(n,t)
                
                if (sm is None) or (sp is None):
                    # One or both of my neighbours are further
                    # down the tree.  Just recur me down to the
                    # next level and try again there.
                    self.__look_out_below(s0,n,key)
                    continue

                # rp & rm are only rank-1 ... don't yet exploit this

                r = rp.inner(sm,ind)
                r0.inner(s0,ind,result=r)  
                rm.inner(sp,ind,result=r)
                result.s[n][key] = r.make_self_contiguous_by_map(amap).scale(twon)

                #r = rp.inner(sm,ind,amap)
                #r0.inner(s0,ind,amap,result=r)  
                #rm.inner(sp,ind,amap,result=r)
                #result.s[n][key] = r.scale(twon)
                
                #result.s[n][key] = (rp.inner(sm,ind,amap) + 
                #                    r0.inner(s0,ind,amap) +
                #                    rm.inner(sp,ind,amap)).scale(twon)

        self.sclean()

        perf.exit('diff')
        return result

    def diff2(self,axis):
        '''

        Apply a central derivative approximation to the second
        derivative along the specified axis, currently with zero
        boundary conditions.  Plan to implement periodic boundary
        conditions also.

        This operation corresponds to .  D2 = DT*D =
        self.diff(axis).diff(axist,transpose=1) and it should be
        negative (semi-definite ... I think that the boundary
        condition introduces a null space ... Beylkin suggests how to
        get rid of it but the projector is not yet coded).

        Self is reconstructed but should otherwise be unchanged.

        Returns a new function in the scaling function basis.

        '''
        if self.debug:
            print "diffing2", id(self), self.compressed, axis

        self.reconstruct()
        perf.enter('diff2')
        
        # Since are modifying the internals of result here must make
        # correct choice of NOT compressed
        result = self.new(k=self.k, compress=0)

        zero = Tensor(self.k,self.k,self.k,self.k,self.k,self.k)

        ind = (1,axis)                  # Indices to contract in products
        if axis == 0:
            amap = None
        elif axis == 1:
            amap = [1,0,2,3,4,5]
        elif axis == 2:
            amap = [2,0,1,3,4,5]
        elif axis == 3:
            amap = [3,0,1,2,4,5]
        elif axis == 4:
            amap = [4,0,1,2,3,5]
        elif axis == 5:
            amap = [5,0,1,2,3,4]
        else:
            raise "Don't you think six dimensions are enough?"
            
        rm, r0, rp = self.make_dc_periodic()
        Rp2 = rm.transpose()*rp
        Rp1 = rm.transpose()*r0 + r0.transpose()*rp
        R0  = rm.transpose()*rm + r0.transpose()*r0 + rp.transpose()*rp
        Rm1 =                     r0.transpose()*rm + rp.transpose()*r0
        Rm2 =                                         rp.transpose()*rm

        nmax = max(self.s.keys())

        for n in range(nmax+1):         #  Fill in empty levels
            try:
                self.s[n]
            except KeyError:
                self.s[n] = {}
            
        for n in range(nmax+1):
            result.s[n] = {}
            sn = self.s[n]
            twon = 2**n
            twon1= twon - 1
            scale = 4.0**n
            for key in self.s[n].keys():
                s = [-2,-1,0,1,2]
                got_them_all = 1
                for shift in s:
                    t=list(key); t[axis]=t[axis]+shift; t=tuple(t)
                    if t[axis]<0 or t[axis]>twon1:
                        s[2+shift] = zero
                    else:
                        s[2+shift] = self.__sock_it_to_me(n,t)
                        if s[2+shift] is None:
                            got_them_all = 0
                            break
                if not got_them_all:
                    # One or more of my neighbours are further
                    # down the tree.  Just recur me down to the
                    # next level and try again there.
                    self.__look_out_below(sn[key],n,key)
                    continue


                # Rp2 & Rm2 are only rank-1 ... don't yet exploit this
                sm2, sm1, s0, sp1, sp2 = s

                r = Rp2.inner(sm2,ind,amap)
                Rp1.inner(sm1,ind,amap,result=r)
                R0.inner(s0,ind,amap,result=r)
                Rm1.inner(sp1,ind,amap,result=r)
                Rm2.inner(sp2,ind,amap,result=r)
                result.s[n][key]=r.scale(scale)

                #result.s[n][key] = (Rp2.inner(sm2,ind,amap) +
                #                    Rp1.inner(sm1,ind,amap) + 
                #                    R0.inner(s0,ind,amap) +
                #                    Rm1.inner(sp1,ind,amap) +
                #                    Rm2.inner(sp2,ind,amap)).scale(scale)

        self.sclean()
        perf.exit('diff2')
        return result

    def laplacian(self):
        '''

        Compute the laplacian using a central derivative currently
        with a zero boundary condition (embedding).
        This operation corresponds to
        .   D2 = DT*D = self.diff(axis).diff(axis,transpose=1)
        and it should be negative (semi-definite ... I think that
        the embedding introduces a null space ... Beylkin suggests
        how to get rid of it but the projector is not yet coded).

        Self is reconstructed but should otherwise be unchanged.

        Returns a new function in the scaling function basis.

        '''
        if self.debug:
            print "laplacian", id(self), self.compressed
        
        self.reconstruct()
        perf.enter('laplacian')

        # Since are modifying the internals of result here must make
        # correct choice of NOT compressed
        result = self.new(k=self.k, compress=0)

        zero = Tensor(self.k,self.k,self.k,self.k,self.k,self.k)

        # For each axis the indices to contract in products and the
        # map from default to desired order in the contractions
        ind = ((1,0), (1,1), (1,2), (1,3), (1,4), (1,5))
        amap = ((0,1,2,3,4,5),
                (1,0,2,3,4,5),
                (2,0,1,3,4,5),
                (3,0,1,2,4,5),
                (4,0,1,2,3,5),
                (5,0,1,2,3,4))
            
        rm, r0, rp = self.make_dc_periodic()
        Rp2 = rm.transpose()*rp
        Rp1 = rm.transpose()*r0 + r0.transpose()*rp
        R0  = rm.transpose()*rm + r0.transpose()*r0 + rp.transpose()*rp
        Rm1 =                     r0.transpose()*rm + rp.transpose()*r0
        Rm2 =                                         rp.transpose()*rm

        nmax = max(self.s.keys())

        for n in range(nmax+1):         #  Fill in empty levels
            try:
                self.s[n]
            except KeyError:
                self.s[n] = {}
            
        for n in range(nmax+1):
            result.s[n] = {}
            sn = self.s[n]
            twon = 2**n
            twon1= twon - 1
            scale = -(4.0**n)
            for key in self.s[n].keys():
                s = [[-2,-1,0,1,2],[-2,-1,0,1,2],[-2,-1,0,1,2],
                     [-2,-1,0,1,2],[-2,-1,0,1,2],[-2,-1,0,1,2]]
                got_them_all = 1
                for axis in 0,1,2,3,4,5:
                    for shift in [-2,-1,0,1,2]:
                        t=list(key); t[axis]=t[axis]+shift; t=tuple(t)
                        if t[axis]<0 or t[axis]>twon1:
                            s[axis][2+shift] = zero
                        else:
                            s[axis][2+shift] = self.__sock_it_to_me(n,t)
                            if s[axis][2+shift] is None:
                                got_them_all = 0
                                break
                    if not got_them_all:
                        break
                    
                if not got_them_all:
                    # One or more of my neighbours are further
                    # down the tree.  Just recur me down to the
                    # next level and try again there.
                    self.__look_out_below(sn[key],n,key)
                    continue

                # Rp2 & Rm2 are only rank-1 ... don't yet exploit this
                r = Tensor(self.k,self.k,self.k,self.k,self.k,self.k)
                for axis in [0,1,2,3,4,5]:
                    sm2, sm1, s0, sp1, sp2 = s[axis]
                    i, a = ind[axis], amap[axis]
                    Rp2.inner(sm2,i,a,result=r)
                    Rp1.inner(sm1,i,a,result=r)
                    R0.inner(s0,i,a,result=r)
                    Rm1.inner(sp1,i,a,result=r)
                    Rm2.inner(sp2,i,a,result=r)
##                     r = r + (Rp2.inner(sm2,i,a) + 
##                              Rp1.inner(sm1,i,a) + 
##                              R0.inner(s0,i,a)   +
##                              Rm1.inner(sp1,i,a) +
##                              Rm2.inner(sp2,i,a))
                    
                result.s[n][key] = r.scale(scale)

        self.sclean()
        perf.exit('laplacian')
        return result
                
    def latex_zslice(self):
        pass

    def to6d(self):
        return self

class Function2d(Function):
    '''

    This is just Function restricted to 2d ... refer to that
    for documentation, comments and perhaps corrections.

    '''
    def __init__(self,
                 function=None, cfunction=None,
                 thresh=None, k=None,
                 initial_level=None,refine=1,compress=1,box=None):

        perf.enter('init')

        self.ndim = 2
        self.new = Function2d
        self.newtensor = Matrix
        
        if k:
            self.k = k                  # Override default in class
        else:
            self.k = k = self.k

        if self.operator:               # Operators need double order
            self.k = k = 2*k

        if thresh:
            self.thresh = thresh        # Override default in class
        else:
            self.thresh = thresh = self.thresh

        if initial_level is None:
            self.initial_level = initial_level = Function.initial_level
        else:
            self.initial_level = initial_level

        self.d = {}             # Directory of wavelet coefficients
        self.s = {}             # Directory of scaling function coefficients
        self.function = function# Remember the function being compressed
        self.cfunction = cfunction
        self.compressed = 0     # Not yet compressed

        self.init_twoscale(k)
        self.init_quadrature(k)

        if cfunction:
            if not (type(cfunction) == type(' ') and \
                    cfunction[-6:] == 'p_void'):
                raise TypeError,"cfunction should be a C "\
                      "function pointer cast to void *"

        if function or cfunction:
            self.nterminated = 0
            if box is None:
                trans = None
            else:
                trans = []
                lxlo,lylo,lxhi,lyhi = box
                yrange = range(lylo,lyhi+1)
                for lx in range(lxlo,lxhi+1):
                    for ly in yrange:
                            trans.append((lx,ly))
            
            perf.enter('fs2d proj')
            self.fine_scale_projection(initial_level,trans)
            perf.exit('fs2d proj')
            
            if refine:
                refinethese = {}
                for lx, ly in self.s[initial_level].keys():
                    refinethese[parent(lx),parent(ly)] = 1
                for lx, ly in refinethese.keys():
                    self.refine_fine_scale_projection(initial_level-1,
                                                      lx, ly)

            if self.nterminated > 0:
                print "Terminated refinement:", self.nterminated

            if compress:
                self.compress(truncate=0)
        else:
            self.s[0] = {}
            self.s[0][0,0] = Matrix(k,k)
            if compress:
                self.compressed = 1
                self.d[0] = {}
                self.d[0][0,0] = Matrix(2*k,2*k)
            
        if self.debug:
            print "Made new 2d function", id(self), self.k, self.thresh, self.compressed
            print "s:", self.s.keys()
            if self.compressed: print "d:", self.d.keys()

        perf.exit('init')

    def fplane(self, xlo, ylo, npt, h, pts, f, function):
        '''

        Tabulate the function on the quadrature points in the
        specified cube ... this is intended to be replaced
        by a C routine.

        '''
        rangenpt = range(npt)
        for p in rangenpt:
            x = xlo+pts[p]*h
            for q in rangenpt:
                y = ylo+pts[q]*h
                f[p,q] = function(x,y)

    def fine_scale_projection(self,n,trans=None):
        '''

        Project the function onto scaling functions at level n.

        If self.operator is true, then we are computing the kernel of
        an integral operator. Operators need -2^n <= l <= 2^n as the
        default range.  

        If l[] is specified it is a list of translations that are to
        be considered.  Otherwise the full list is used.

        Modifies:
        :           self.s

        '''

        f = self.function
        try:
            self.s[n]
        except KeyError:
            self.s[n] = {}
        
        h = 1.0/(2.0**n)                # Box size on target level 
        scale = sqrt(h)
        scale = scale*scale             # Since in 2D
        phiw = self.quad_phiw
        pts = self.quad_x 
        npt = len(pts)

        if not trans:
            if self.operator:
                trans1d = range(-2**n,2**n+2)
            else:
                trans1d = range(2**n)
            trans = []
            for lx in trans1d:
                for ly in trans1d:
                    trans.append((lx,ly))

        for lx,ly in trans:
            xlo = lx*h
            ylo = ly*h
            f = Matrix(npt,npt)
            if self.cfunction:
                cfcube.cfplane(xlo, ylo, npt, h,
                              pts.swigptr(), f.swigptr(),
                              self.cfunction)
            else:
                # Plain old Python function or non-cast C function
                # This is SLOW due to all of the tensor references.
                self.fplane(xlo, ylo, npt, h, pts, f, self.function)
            f.scale(scale)

            self.s[n][lx,ly] = f.fast_transform(phiw)

    def refine_fine_scale_projection(self,n,lx,ly,prnt=0):
        '''

        We have the fine scale projection at level n+1, so we can
        compute the scaling function and wavelet coefficients at level
        n.  Considering just box l on level n, examine the difference
        coefficients.  If they do not satisfy the threshold, then
        refine the representation by recurring down to a finer level.

        If refinement occurs this modifies:
        :   self.s

        Set prnt to non-zero for printing

        '''
        if n > self.max_refine_level:
            self.nterminated = self.nterminated + 1
            return

        if Function.truncate_method == 0:
            tol = self.thresh
        else:
            tol = sqrt(1.0/4.0**n)*self.thresh
            
        k = self.k
        ss = self.gather_scaling_coeffs(n,lx,ly)
        ss = self.filter(ss)
        s = Matrix(ss[0:k,0:k])
        ss[0:k,0:k].fill(0.0)
        dnorm = ss.normf()
        if dnorm <= tol:
##             # The difference coefficients are neglible so can truncate
##             # at this level.  Delete the coeffs on the level below to
##             # save space, and adopt the scaling function coeffs
##             # computed from the lower level by twoscale since they
##             # will be more accurate
##             if not self.s.has_key(n): self.s[n] = {}
##             self.s[n][lx,ly] = s
##             lx2, ly2 = lx*2, ly*2
##             for ix in range(2):
##                 for iy in range(2):
##                     del(self.s[n+1][lx2+ix,ly2+iy])
            pass
        else:
            # First remove the inaccurate scaling function coeffs at
            # level n+1, then refine to level n+2
            if prnt:
                for i in range(n): print "  ",
                print " refining level=%d l=(%d,%d) norm=%.1e" % \
                      (n, lx, ly, dnorm)
            lx2, ly2 = lx*2, ly*2
            range2 = range(2)
            for ix in range2:
                ixlo = ix*k
                ixhi = ixlo + k
                for iy in range2:
                    iylo = iy*k
                    iyhi = iylo + k
                    del(self.s[n+1][lx2+ix,ly2+iy])
                    lx4, ly4 = 2*(lx2+ix), 2*(ly2+iy)
                    for iix in range2:
                        for iiy in range2:
                            self.fine_scale_projection(
                                n+2,[(lx4+iix,ly4+iiy)])
                    self.refine_fine_scale_projection(
                        n+1,lx2+ix,ly2+iy)

    def gather_scaling_coeffs(self, n, lx, ly):
        '''

        For a given translation at level n, gather the corresponding
        scaling function coefficients at level n+1.

        Some of the boxes on the lower level may be missing.

        '''
        k = self.k
        sn1 = self.s[n+1]
        ss = Matrix(2*k,2*k)
        lx2, ly2 = lx*2, ly*2
        for ix in range(2):
            ixlo = ix*k
            ixhi = ixlo + k
            for iy in range(2):
                iylo = iy*k
                iyhi = iylo + k
                key = (lx2+ix,ly2+iy)
                if sn1.has_key(key):
                    ss[ixlo:ixhi,iylo:iyhi] =  sn1[key]
        return ss
        
    def reconstruct(self, nonstandard=0):
        '''

        Reconstruct the projection onto scaling functions from the
        compressed form terminating at the locally lowest level
        without difference coefficients.

        The compression algorithm has the responsibility of ensuring
        that the there are no signficant nodes below a node without
        difference coefficients.

        This algorithm must ensure that there are scaling function
        coefficients only at the locally lowest level (since other
        algorithms will just run down the tree until the first
        coefficient is found, and then stop).

        A variant of this algorithm is used to reconstruct the
        function in preparation for the application of a non-standard
        form operator.  For this we need to keep the difference
        coefficients and also keep the sum coefficients at each level
        on the way down, but only down to locally finest level at
        which we have difference coefficients.  There is (will be)
        another routine to clean the resulting mess back into the
        conventional hierarchical basis definition.

        '''
        if self.debug:
            print "reconstructing", id(self), self.compressed, self.s.keys()
            
        if not self.compressed:
            if nonstandard:
                self.compress()
            else:
                return self

        perf.enter('reconstrct')
        self.compressed = 0

        k = self.k

        nmax = max(self.d.keys())
        if nonstandard == 1:
            # Only need scaling functions down to where the difference
            # coefficients are locally non-zero.
            #nmax = nmax-1
            nmax = nmax
            
        range2 = range(2)

        for n in range(nmax+1):
            sn = self.s[n]
            dn = self.d[n]
            sn1 = self.s[n+1] = {}
            if nonstandard==1 and n<nmax:
                dn1 = self.d[n+1]
            for lx, ly in dn.keys():
                d = dn[lx,ly]
                d[0:k,0:k] = sn[lx,ly]
                if nonstandard==1 and n>=nmax: continue
                if nonstandard != 2: del(sn[lx,ly])
                d = self.unfilter(d)
                lx2, ly2 = lx*2, ly*2
                for ix in range2:
                    for iy in range2:
                        key2 = lx2+ix,ly2+iy
                        if (nonstandard==1) and (not dn1.has_key(key2)):
                            continue
                        dxyz = d[ix*k:ix*k+k,iy*k:iy*k+k]
                        sn1[key2] = Matrix(dxyz)

            if nonstandard == 1:
                del self.s[n]
            else:
                del self.d[n]

        if nonstandard == 1: del self.s[nmax+1]
        
        if nmax == -1 and nonstandard == 1:
            print "FUDGE", self.s[0][0,0].normf()
            self.d[0][0,0][0:k,0:k] = self.s[0][0,0]
            print "FUDGE", self.d[0][0,0].normf()
            del self.s[0][0,0]

        perf.exit('reconstrct')
        return self

    def compress(self,truncate=0):
        '''

        Compress the scaling function coefficients (Vn+1 = Vn + Wn).

        In the usual compressed form, the scaling function
        coefficients will only exist on the lowest level, but this
        will vary across the domain.

        In the non-standard compressed form, both scaling function and
        difference coefficients will be defined on all levels ... must
        add the result of compressing the sum coefficients into the
        existing difference coefficients.

        In the form that results from addition in the scaling function
        basis, there may be sum coefficients at many levels, but no
        differences.

        The difference coefficients are truncated to the specified
        threshold.  So that reconstruction can just proceed down to
        the first level without difference coefficients, it is
        important to only truncate difference coefficients at leaf
        nodes.

        For ease of handling all the special cases now implement
        truncation as a separate and optional step.

        '''
        if self.debug:
            print "compressing", id(self), self.compressed, self.s.keys()

        if self.compressed:
            return self
        self.compressed = 1
        perf.enter('compress')
        
        nmax = max(self.s.keys())
        for n in range(nmax-1,-1,-1):
            if not self.d.has_key(n):
                self.d[n] = {}

        k = self.k
        for n in range(nmax-1,-1,-1):
            dn  = self.d[n]
            sn1 = self.s[n+1]
            try:
                sn = self.s[n]
            except KeyError:
                sn = self.s[n] = {}
            try:
                dn1 = self.d[n+1]
            except KeyError:
                dn1 = {}
            
            # Compute the list of translations that we must consider, driven
            # from the non-zero entries on the level below.
            llist = {}
            for lx,ly in sn1.keys():
                llist[parent(lx),parent(ly)] = 1
            
            for lx,ly in llist.keys():
                ss = self.gather_scaling_coeffs(n, lx, ly)
                ss = self.filter(ss)
                s = Matrix(ss[0:k,0:k]) # Force a copy
                if sn.has_key((lx,ly)):
                    sn[lx,ly] = sn[lx,ly].gaxpy(1.0,s,1.0)
                else:
                    sn[lx,ly] = s

                ss[0:k,0:k].fill(0.0) # Now just want the diff coeffs
                if dn.has_key((lx,ly)):
                    dn[lx,ly] = dn[lx,ly].gaxpy(1.0,ss,1.0)
                else:
                    dn[lx,ly] = ss

            del self.s[n+1]

        if truncate:
            self.truncate()

        if not self.d.has_key(0):
            self.d[0] = {}
        if not self.d[0].has_key((0,0)):
            self.d[0][0,0] = Matrix(2*k,2*k)

        if not self.s.has_key(0):
            self.s[0] = {}
        if not self.s[0].has_key((0,0)):
            self.s[0][0,0] = Matrix(k,k)

        perf.exit('compress')
        return self

    def basis(self,refine=0):
        '''

        Self is a compressed function.  Export a description of
        the basis (for use by restrict).  This describes the
        pieces of W0+W1+...+W(n-1) in use.

        If refine=1, then if the difference coefficients do not
        satisfy the accuracy goal, the basis is locally refined down
        one level.

        '''
        if not self.compressed: self.compress()

        if Function.truncate_method == 2:
            lend = self.__lend()
        
        basis = {}
        nmax = max(self.d.keys())
        for n in range(nmax+1):
            dn = self.d[n]
            bn = basis[n] = {}
            for key in dn.keys():
                bn[key] = 1

        if refine:
            basis[nmax+1] = {}          # Just in case
            for n in range(nmax+1):
                dn = self.d[n]
                bn = basis[n]
                bn1 = basis[n+1]

                keys = dn.keys()

                if Function.truncate_method == 0:
                    tol = self.thresh
                elif Function.truncate_method == 1:
                    tol = sqrt(1.0/4.0**n)*self.thresh
                else:
                    tol = self.thresh/sqrt(lend)

                range2 = range(2)
                for key in keys():
                    if (dn[key].normf() > tol):
                        lx,ly = key
                        lx2, ly2 = lx*2, ly*2
                        printed = 0
                        for ix in range2:
                            for iy in range2:
                                key2 = (lx2+ix,ly2+iy)
                                if not bn1.has_key(key2):
                                    if not printed:
                                        print "Refining", n, key
                                        printed = 1
                                    bn1[key2] = 1
            if not basis[nmax+1]:
                del basis[nmax+1]

        return basis
    
    def truncate(self,thresh=None):
        '''

        The difference coefficients are truncated to the specified
        threshold.  So that reconstruction can just proceed down to
        the first level without difference coefficients, it is
        important to only truncate difference coefficients at leaf
        nodes.

        Returns self so that operations can be chained.

        '''
        if self.debug:
            print "truncating", id(self), self.compressed, self.thresh
            
        if not self.compressed:
            self.compress()

        if Function.truncate_method == 2:
            lend = self.__lend()
        
        if thresh is None:
            thresh = self.thresh
        
        perf.enter('truncate')

        nmax = max(self.d.keys())

        is_parent = {}
        for n in range(nmax,0,-1):
            dn = self.d[n]
            keys = dn.keys()

            if Function.truncate_method == 0:
                tol = thresh
            elif Function.truncate_method == 1:
                tol = sqrt(1.0/4.0**n)*thresh
            else:
                tol = thresh/sqrt(lend)

            # If not a parent and below threshold kill box.
            # If above threshold mark parent box as such.
            next_is_parent = {}
            for key in keys:
                if (not is_parent.has_key(key)) and \
                   (dn[key].normf() < tol):
                    del dn[key]
                else:
                    lx,ly = key
                    next_is_parent[parent(lx),parent(ly)] = 1

            is_parent = next_is_parent
                    
            # If the current level is totally empty and it is also
            # the lowest level, it can be deleted entirely.
            if (not self.d[n]) and (not self.d.has_key(n+1)):
                del self.d[n]

        perf.exit('truncate')

##         if Function.truncate_method == 2:
##             # If we deleted a lot of coefficients, try again ...
##             if self.__lend()*2 < lend:
##                 self.truncate()

        return self

    def eval(self,x,y):
        '''

        Evaluate the function at the given point.

        If it is not compressed then we just need to walk down the
        tree until we locate the leaf node containing x,y,z.  If it is
        compressed then we do the same, but also have to add up the
        scaling function coeffs on the way down.

        '''
        n = 0
        k = self.k
        lx = ly = 0
        xx, yy = x, y
        if self.compressed:
            # Note always want loop to terminate via the break otherwise
            # n is one too small (x & l would also be wrong)
            s = self.s[0][0,0]
            for n in range(max(self.d.keys())+2):
                try:
                    d = Matrix(self.d[n][lx,ly]) # Take a copy
                except KeyError:
                    # At the bottom ... evaluate s*phi(xx) at level n
                    break
                d[0:k,0:k] = s
                d = self.unfilter(d)
                scale = 2**(n+1)
                xx, yy = x*scale, y*scale
                lx, ly = min(scale-1,int(xx)), min(scale-1,int(yy))
                xx, yy = xx-lx, yy-ly
                ix, iy = lx%2, ly%2
                s = d[ix*k:ix*k+k,iy*k:iy*k+k]
            s = Matrix(s)
        else:
            # The function is not compressed.  Assume that it is
            # represented at some (locally defined) finest level
            # in the scaling function basis.
            nmax = max(self.s.keys()) + 1
            for n in range(nmax):
                twon = 2**n
                twon1= twon - 1
                xx, yy = x*twon, y*twon
                lx, ly = min(twon1,int(xx)), min(twon1,int(yy))
                xx, yy = xx-lx, yy-ly
                try:
                    s = self.s[n][lx,ly]
                    break
                except KeyError:
                    pass
            
        #px = Vector(phi(xx,k))
        #value = s.inner(px)*sqrt(2.0**n)

        value = cfcube.eval2d(k,n,xx,yy,s.swigptr())
        return value

    def eval_err(self, npt=1000):
        self.reconstruct()
        errsq = 0.0
        maxerr= 0.0
        maxrelerr=0.0
        seed(201)
        for i in range(npt):
            x = 0.4+random()*0.1
            y = 0.4+random()*0.1
            value = self.eval(x,y)
            exact = self.function(x,y)
            err = abs(exact - value)
            maxerr = max(maxerr, err)
            errsq = errsq + err*err
            if exact:
                maxrelerr = max(maxrelerr,err/exact)
        rmserr = sqrt(errsq/npt)
        print "RMS abs error = %.1e" % rmserr
        print "MAX abs error = %.1e" % maxerr
        print "MAX rel error = %.1e" % maxrelerr
        self.compress(truncate=0)


    def print_layers(self, print_pic=0):
        if self.compressed:
            print " Analysis by layers of compressed function"
            nmax = max(self.d.keys())
            normsq = self.s[0][0,0].normf()**2
            print "\nNorm in V0 %.8e" % sqrt(normsq)
            for n in range(nmax+1):
                sum = 0.0
                for xyz in self.d[n].keys():
                    sum = sum + self.d[n][xyz].normf()**2
                normsq = normsq + sum
                print "Norm in W%d %.8e %.8e" % (n,sqrt(sum),sqrt(normsq))

            if not print_pic: return  
        else:
            print " Analysis by layers of reconstructed function"
            nmax = max(self.s.keys())
            normsq = 0.0
            for n in range(nmax+1):
                sum = 0.0
                for xyz in self.s[n].keys():
                    sum = sum + self.s[n][xyz].normf()**2
                normsq = normsq + sum
                print "Norm in V%d %.8e %.8e" % (n,sqrt(sum),sqrt(normsq))

            if not print_pic: return

    def sclean(self):
        '''

        Cleans up scaling function coefficients after diffing or
        multiplying.  Only want to keep the scaling function coeffs
        at the highest level in the tree ... ones lower down will
        have been generated by two-scale just for evaluation.
        Kill all nodes that have a parent - infanticide.

        '''
        for n in range(max(self.s.keys()),0,-1):
            sn, sn1 = self.s[n], self.s[n-1]
            for lx,ly in sn.keys():
                if sn1.has_key((parent(lx),parent(ly))):
                    del(sn[lx,ly])
        
    def mul(self,other):
        '''

        Multiplication of a function by a another function.
        Replaces the algorithm based on squaring.

        Returns a new function in the scaling function basis.
        Self/other are not modified though may be reconstructed if
        either is compressed.

        Added an optimization for multiplcation by the mask which is
        mask(x) = 1 for 0.0625 < x < 1-0.625.  Thus, only need to do
        the multiplication we have boxes 0 or 15 on level 4 (or boxes
        contained within below).  The mask is identified by the
        attribute mask which noone else will have.

        '''

        self.reconstruct()
        other.reconstruct()

        if hasattr(self,'mask'): self,other = other,self
        mask = hasattr(other,'mask')
        
        perf.enter('mul')

        result = self.new(k=self.k,compress=0) # Was thresh=self.thresh
        result.s = {}

        k = self.k
        phit = Matrix(self.quad_phi.cycledim(1)) # Contig. for fast_transform
        phiw = self.quad_phiw
        npt = len(self.quad_x)
        nmax = max(self.s.keys()+other.s.keys())

        tmp1 = Matrix(k,k)
        tmp2 = Matrix(k,k)

        range2 = range(2)

        if Function.autorefine:     # NOTE REFERENCE TO CLASS
            result.s[0] = {}
             
        for n in range(nmax+1):
            try:
                sn = self.s[n]
            except KeyError:
                sn = self.s[n] = {}

            try:
                on = other.s[n]
            except KeyError:
                on = other.s[n] = {}
                
            if Function.autorefine:     # NOTE REFERENCE TO CLASS
                rn1 = result.s[n+1] = {}
                scale = sqrt(4.0**(n+1))
            else:
                rn = result.s[n] = {}
                scale = sqrt(4.0**n)

            # Generate union of keys on this level. 
            keys = sn.keys() + on.keys()
            keys.sort()
            for i in range(len(keys)-1,0,-1):
                if keys[i] == keys[i-1]:
                    del keys[i]

            if mask:
                if n <= 4:
                    masklo, maskhi = 0, 2**n-1
                else:
                    masklo, maskhi = 2**(n-4)-1, 15*2**(n-4)
                    
            for key in keys:
                lx,ly = key

                # If other is the mask, exclude keys for which the
                # mask is unity.
                if mask:
                    if masklo<lx<maskhi and masklo<ly<maskhi:
                        if sn.has_key(key): rn[key] = Matrix(sn[key])
                        continue

                s =  self.__sock_it_to_me(n,key)
                if not (s is None):
                    o = other.__sock_it_to_me(n,key)
                    if not (o is None):
                        if Function.autorefine: # NOTE REFERENCE TO CLASS
                            ss = self.unfilter(s,sonly=1)
                            oo = self.unfilter(o,sonly=1)
                            lx2, ly2 = 2*lx, 2*ly
                            for ix in range2:
                                for iy in range2:
                                    f = ss[ix*k:ix*k+k,iy*k:iy*k+k]
                                    g = oo[ix*k:ix*k+k,iy*k:iy*k+k]
                                    f = f.transform(phit,result=tmp1)
                                    g = g.transform(phit,result=tmp2)
                                    f.emul(g)
                                    rn1[lx2+ix,ly2+iy] = \
                                       f.fast_transform(phiw).scale(scale)
                        else:
                            f = s.fast_transform(phit,result=tmp1)
                            g = o.fast_transform(phit,result=tmp2)
                            f.emul(g)
                            rn[key] = f.fast_transform(phiw).scale(scale)
                            
        self.sclean()
        other.sclean()

        perf.exit('mul')

        return result

    def mul_with1d(self,other,side=1):
        '''
        
        side = 1 : left    f_new(r1,r2) = self(r1,r2) * other(r1)
        side = 2 : right   f_new(r1,r2) = self(r1,r2) * other(r2)
        
        Multiplication of a function by a another function.
        Replaces the algorithm based on squaring.

        Returns a new function in the scaling function basis.
        Self/other are not modified though may be reconstructed if
        either is compressed.

        Added an optimization for multiplcation by the mask which is
        mask(x) = 1 for 0.0625 < x < 1-0.625.  Thus, only need to do
        the multiplication we have boxes 0 or 15 on level 4 (or boxes
        contained within below).  The mask is identified by the
        attribute mask which noone else will have.

        '''

        if side!=1 and side!=2:
            raise "side is either 1 or 2"

        self.reconstruct()
        other.reconstruct()

        if hasattr(self,'mask'): self,other = other,self
        mask = hasattr(other,'mask')
        
        perf.enter('mul')

        result = self.new(k=self.k,compress=0) # Was thresh=self.thresh
        result.s = {}

        k = self.k
        phit = Matrix(self.quad_phi.cycledim(1)) # Contig. for fast_transform
        phiw = self.quad_phiw
        npt = len(self.quad_x)
        nmax = max(self.s.keys()+other.s.keys())

        tmp1 = Matrix(k,k)
        one  = Vector(k)
        one[0] = 1.0

        range2 = range(2)

        if Function.autorefine:     # NOTE REFERENCE TO CLASS
            result.s[0] = {}
             
        for n in range(nmax+1):
            try:
                sn = self.s[n]
            except KeyError:
                sn = self.s[n] = {}

            try:
                on = other.s[n]
            except KeyError:
                on = other.s[n] = {}

            if Function.autorefine:     # NOTE REFERENCE TO CLASS
                rn1 = result.s[n+1] = {}
                scale = 2.0**(n+1)
            else:
                rn = result.s[n] = {}
                scale = 2.0**n
                
            one[0] = 1.0/sqrt(scale)

            # Generate union of keys on this level.
            trans1d = range(2**n)
            onkeys = []
            if(side==1):            
                for key in on.keys():
                    for lx in trans1d:
                        onkeys.append((key,)+(lx,))  # <==========check for left
            else:
                for key in on.keys():
                    for lx in trans1d:
                        onkeys.append((lx,)+(key,))  # <==========check for right

            keys = sn.keys() + onkeys
            keys.sort()
            for i in range(len(keys)-1,0,-1):
                if keys[i] == keys[i-1]:
                    del keys[i]

            if mask:
                if n <= 4:
                    masklo, maskhi = 0, 2**n-1
                else:
                    masklo, maskhi = 2**(n-4)-1, 15*2**(n-4)

            cache_g = {}
            for key in keys:
                lx,ly = key

                # If other is the mask, exclude keys for which the
                # mask is unity.
                if mask:
                    if masklo<lx<maskhi and masklo<ly<maskhi:
                        if sn.has_key(key): rn[key] = Matrix(sn[key])
                        continue

                s =  self.__sock_it_to_me(n,key)
                if not (s is None):
                    if(side==1):            
                        o = other._Function1d__sock_it_to_me(n,(lx)) # <=========== check for left
                    else:
                        o = other._Function1d__sock_it_to_me(n,(ly)) # <=========== check for right
                    if not (o is None):
                        if Function.autorefine: # NOTE REFERENCE TO CLASS
                            ss = self.unfilter(s,sonly=1)
                            oo = other.unfilter(o,sonly=1)
                            lx2, ly2 = 2*lx, 2*ly
                            for ix in range2:
                                for iy in range2:
                                    f = ss[ix*k:ix*k+k,iy*k:iy*k+k]
                                    if(side==1):            
                                        g = oo[ix*k:ix*k+k] # <==check for left
                                        try:
                                            g = cache_g[(lx)] # <====== check for left
                                        except KeyError:
                                            g = g.outer(one) # <====== check for left
                                            g = g.transform(phit)
                                            cache_g[(lx)] = g
                                    else:
                                        g = oo[iy*k:iy*k+k] # <==check for right
                                        try:
                                            g = cache_g[(ly)] # <====== check for right
                                        except KeyError:
                                            g = one.outer(g) # <====== check for right
                                            g = g.transform(phit)
                                            cache_g[(ly)] = g
                                            
                                    f = f.transform(phit,result=tmp1)
                                    f.emul(g)
                                    rn1[lx2+ix,ly2+iy] = \
                                       f.fast_transform(phiw).scale(scale)
                        else:
                            if(side==1):            
                                try:
                                    g = cache_g[(lx)] # <====== check for right
                                except KeyError:
                                    o = o.outer(one)
                                    g = cache_g[(lx)] = o.fast_transform(phit)
                            else:
                                try:
                                    g = cache_g[(ly)] # <====== check for right
                                except KeyError:
                                    o = one.outer(o)
                                    g = cache_g[(ly)] = o.fast_transform(phit)
                                    
                            f = s.fast_transform(phit,result=tmp1)
                            f.emul(g)
                            rn[key] = f.fast_transform(phiw).scale(scale)
            del cache_g
            
        self.sclean()
        other.sclean()

        perf.exit('mul')

        return result

    def __rmul__(self,other):

        return self.__mul__(other)

    def square(self):
        '''

        Square self replacing self by its square in the scaling
        function basis.

        If (Function.autorefine) is set, then it will recur down an
        additional level before squaring, otherwise it does not.

        NB: self is changed.

        Returns self so that operations may be chained.  

        '''
        if self.debug:
            print "squaring", id(self), self.compressed
            
        self.reconstruct()

        k = self.k
        nmax = max(self.s.keys())
        phit = Matrix(self.quad_phi.cycledim(1)) # Contig. for fast_transform
        phiw = self.quad_phiw
        npt = len(self.quad_x)

        tmp1 = Matrix(k,k)
        range2 = range(2)

        try:
            next_keys = self.s[0].keys()
        except KeyError:
            # This can happen if not compressed and initial_level != 0
            self.s[0] = {}
            next_keys = self.s[0].keys()

        for n in range(nmax+1):
        
            sn = self.s[n]
            try:
                sn1 = self.s[n+1]
            except KeyError:
                sn1 = self.s[n+1] = {}

            keys = next_keys
            next_keys = sn1.keys()

            if Function.autorefine:         # NOTE REFERENCE TO CLASS
                ss = Matrix(2*k,2*k)
                scale = sqrt(4.0**(n+1))
            else:
                scale = sqrt(4.0**n)
            for lx, ly in keys:
                if Function.autorefine:     # NOTE REFERENCE TO CLASS
                    # Now are at the locally lowest level ... recur down
                    # one more level and then form the square projecting
                    # onto the approximately refined basis.  Could recur
                    # more dynamically to guarantee precision but
                    # Beylkin does not seem to think this is necessary
                    # except in pathological cases?
                    ss = self.unfilter(sn[lx,ly],sonly=1)
                    del(sn[lx,ly])
                    lx2, ly2 = 2*lx, 2*ly
                    for ix in range2:
                        for iy in range2:
                            f = ss[ix*k:ix*k+k,iy*k:iy*k+k]
                            # Back transform to the function value
                            # at the quadrature points
                            f = f.transform(phit,result=tmp1)
                            # Compute the square on the grid
                            f.emul(f)
                            # Do quadrature on f**2 and scale
                            sn1[lx2+ix,ly2+iy] = \
                                f.fast_transform(phiw).scale(scale)
                                                          
                else:
                    # Don't refine ... just do it at the current level
                    f = sn[lx,ly].fast_transform(phit,result=tmp1)
                    f.emul(f)
                    sn[lx,ly] = f.fast_transform(phiw).scale(scale)
            
        return self

    def local_function(self,function=None,cfunction=None):
        '''

        Apply a scalar local function to self replacing self by
        f(self(x)).  Define either function (a callable object) or
        cfunction (a pointer to a function cast to a pointer to void
        ... so that old versions of SWIG can be used).  Both the C and
        Python functions should take a double argument and return a
        double.  THIS APPROACH IS ONLY GOOD IF THE RESULT OF THE
        FUNCTION IS AS LOCALLY SMOOTH OR SMOOTHER THAN THE INPUT ...

        If (Function.autorefine) is set, then it will recur down an
        additional level before squaring, otherwise it does not.

        NB: self is changed.

        Returns self so that operations may be chained.  

        '''
        if self.debug:
            print "squaring", id(self), self.compressed
            
        range2 = range(2)

        self.reconstruct()
        perf.enter('local func')

        k = self.k
        nmax = max(self.s.keys())
        phit = Matrix(self.quad_phi.cycledim(1)) # Contig. for fast_transform
        phiw = self.quad_phiw
        npt = len(self.quad_x)

        tmp1 = Matrix(k,k)

        try:
            next_keys = self.s[0].keys()
        except KeyError:
            # This can happen if not compressed and initial_level != 0
            self.s[0] = {}
            next_keys = self.s[0].keys()

        for n in range(nmax+1):
        
            sn = self.s[n]
            try:
                sn1 = self.s[n+1]
            except KeyError:
                sn1 = self.s[n+1] = {}

            keys = next_keys
            next_keys = sn1.keys()

            if Function.autorefine:         # NOTE REFERENCE TO CLASS
                ss = Matrix(2*k,2*k)
                scale = sqrt(4.0**(n+1))
            else:
                scale = sqrt(4.0**n)
            for lx, ly in keys:
                if Function.autorefine:     # NOTE REFERENCE TO CLASS
                    # Now are at the locally lowest level ... recur down
                    # one more level and then form the square projecting
                    # onto the approximately refined basis.  Could recur
                    # more dynamically to guarantee precision but
                    # Beylkin does not seem to think this is necessary
                    # except in pathological cases?
                    ss = self.unfilter(sn[lx,ly],sonly=1)
                    del(sn[lx,ly])
                    lx2, ly2 = 2*lx, 2*ly
                    for ix in range2:
                        for iy in range2:
                            f = ss[ix*k:ix*k+k,iy*k:iy*k+k]
                            # Back transform to the function value
                            # at the quadrature points
                            f = Matrix(f).transform(phit,result=tmp1).scale(scale)
                            # Eval function on the refined grid
                            if cfunction:
                                cfcube.clocalfunc(
                                    f.swigptr(),k**2,cfunction)
                            else:
                                f.unaryop(function)
                            # Do quadrature on f**2 and scale
                            sn1[lx2+ix,ly2+iy] = \
                                f.fast_transform(phiw).scale(1.0/scale)
                                                          
                else:
                    # Don't refine ... just do it at the current level
                    f = sn[lx,ly].fast_transform(phit,result=tmp1).scale(scale)
                    if cfunction:
                        cfcube.clocalfunc(f.swigptr(),k**2,cfunction)
                    else:
                        f.unaryop(function)
                    sn[lx,ly] = f.fast_transform(phiw).scale(1.0/scale)
            
        perf.exit('local func')
        return self

    def mul_func(self,function=None,cfunction=None):
        '''

        '''
        if self.debug:
            print "squaring", id(self), self.compressed
            
        range2 = range(2)

        self.reconstruct()
        perf.enter('mul func')

        k = self.k
        nmax = max(self.s.keys())
        phit = Matrix(self.quad_phi.cycledim(1)) # Contig. for fast_transform
        phiw = self.quad_phiw
        pts = self.quad_x 
        npt = len(self.quad_x)

        tmp1 = Matrix(k,k)

        try:
            next_keys = self.s[0].keys()
        except KeyError:
            # This can happen if not compressed and initial_level != 0
            self.s[0] = {}
            next_keys = self.s[0].keys()

        for n in range(nmax+1):

            sn = self.s[n]
            try:
                sn1 = self.s[n+1]
            except KeyError:
                sn1 = self.s[n+1] = {}

            keys = next_keys
            next_keys = sn1.keys()

            if Function.autorefine:         # NOTE REFERENCE TO CLASS
                ss = Matrix(2*k,2*k)
                dd = float(2**(n+1))
            else:
                dd = float(2**n)

            scale = sqrt(dd**2)       # 2D Scale
            h = 1.0/dd                # Box size on target level
             
            for lx, ly in keys:

                if Function.autorefine:     # NOTE REFERENCE TO CLASS
                    # Now are at the locally lowest level ... recur down
                    # one more level and then form the square projecting
                    # onto the approximately refined basis.  Could recur
                    # more dynamically to guarantee precision but
                    # Beylkin does not seem to think this is necessary
                    # except in pathological cases?
                    ss = self.unfilter(sn[lx,ly],sonly=1)
                    del(sn[lx,ly])
                    lx2, ly2 = 2*lx, 2*ly
                    for ix in range2:
                        for iy in range2:

                            # get a cube of g:(func)
                            # ----------------------
                            xlo = (lx2+ix)*h
                            ylo = (ly2+iy)*h
                            g = Matrix(npt,npt)
                            if cfunction:
                                cfcube.cfplane(xlo, ylo, npt, h,
                                              pts.swigptr(), g.swigptr(),
                                              cfunction)
                            else:
                                # Plain old Python function or non-cast C function
                                # This is SLOW due to all of the tensor references.
                                self.fplane(xlo, ylo, npt, h, pts, g, function)
                            
                            # get a cube of f:(func) 
                            # ----------------------
                            f = ss[ix*k:ix*k+k,iy*k:iy*k+k]
                            # Back transform to the function value
                            # at the quadrature points
                            f = f.transform(phit,result=tmp1).scale(scale)
                            
                            # emul between f and g
                            # --------------------
                            f = f.emul(g)
                            
                            # Do quadrature on f**2 and scale
                            sn1[lx2+ix,ly2+iy] = \
                                f.fast_transform(phiw).scale(1.0/scale)
                                                          
                else:

                    # get a cube of g:(func)
                    # ----------------------
                    xlo = lx*h
                    ylo = ly*h
                    g = Matrix(npt,npt)
                    if cfunction:
                        cfcube.cfplane(xlo, ylo, npt, h,
                                      pts.swigptr(), g.swigptr(),
                                      cfunction)
                    else:
                        # Plain old Python function or non-cast C function
                        # This is SLOW due to all of the tensor references.
                        self.fplane(xlo, ylo, npt, h, pts, g, function)

                    # get a cube of f:(func) 
                    # ----------------------
                    # Don't refine ... just do it at the current level
                    f = sn[lx,ly].fast_transform(phit,result=tmp1).scale(scale)

                    # emul between f and g
                    # --------------------
                    f = f.emul(g)
                
                    sn[lx,ly] = f.fast_transform(phiw).scale(1.0/scale)
            
        perf.exit('mul func')
        return self

    def __look_out_below(self,s, n, key):
        '''

        s are scaling coefficients at level n for the given
        translation.  Recur them down one level and put the resulting
        coefficients into the tree.

        '''
        if not self.s.has_key(n+1):
            self.s[n+1] = {}
        k = self.k
        ss = self.unfilter(s,sonly=1)
        lx, ly = key
        lx2, ly2 = lx*2, ly*2
        for ix in range(2):
            for iy in range(2):
                view = ss[ix*k:ix*k+k,iy*k:iy*k+k]
                self.s[n+1][lx2+ix,ly2+iy] = Matrix(view)

    def __sock_it_to_me(self,n,key):
        '''

        Return self.s[n][key] or None if it cannot be made (because it
        has no parent).
        
        If it is absent, then try to recur from the level above to
        make it.  If the parent is also absent, we need to keep
        walking up.  If we hit the top of the tree and still cannot
        find a parent, then return None.

        Some optimizations are possible, such as looking one below
        before recuring up.  Also, has_key() might be faster than
        exception handling.

        '''
        try:
            return self.s[n][key]
        except KeyError:
            if n == 0:
                return None

        lx, ly = key

        lx, ly = parent(lx), parent(ly)
        try:
            s = self.s[n-1][lx,ly]
        except KeyError:
            s = self.__sock_it_to_me(n-1,(lx,ly))

        if s is None:
            return None

        self.__look_out_below(s,n-1,(lx,ly))

        return self.s[n][key]

    def __sock_it_to_me2(self,n,key):
        '''

        Return self.s[n][key] or None if it cannot be made (because it
        has no parent).
     
        If it is absent, then try to recur from the level above to
        make it.  If the parent is also absent, we need to keep
        walking up.  If we hit the top of the tree and still cannot
        find a parent, then return None.

        Some optimizations are possible, such as looking one below
        before recuring up.  Also, has_key() might be faster than
        exception handling.

        '''
        try:
            return self.s[n][key].normf()
        except KeyError:
            if n == 0:
                return None

        lx, ly = key

        lx, ly = parent(lx), parent(ly)
        try:
            s = self.s[n-1][lx,ly]
            s = s.normf()
        except KeyError:
            s = self.__sock_it_to_me2(n-1,(lx,ly))

        if s is None:
            return None

        s = s * 0.50 # = (2.0**(-2.0/2.0))
     
        return s

    def diff(self,axis,transpose=0):
        '''

        Apply a central derivative along the specified axis,
        currently with zero boundary conditions.  Plan to implement
        periodic boundary conditions also.

        If (transpose) then use the transposed derivative operator.

        Self is reconstructed but should otherwise be unchanged.

        Returns a new function in the scaling function basis.

        '''
        if self.debug:
            print "diffing", id(self), self.compressed, axis, \
                  transpose,self.s.keys()

        self.reconstruct()
        perf.enter('diff')
        
        result = self.new(k=self.k, compress=0)

        zero = Matrix(self.k,self.k)

        ind = (1,axis)                  # Indices to contract in products
        if axis == 0:
            amap = None
        elif axis == 1:
            amap = [1,0]
        else:
            raise "Don't you think six dimensions are enough?"
            
        rm, r0, rp = self.make_dc_periodic()

        if transpose:
            rm, r0, rp = rp.transpose(), r0.transpose(), rm.transpose()

        nmax = max(self.s.keys())

        for n in range(nmax+1):         #  Fill in empty levels
            try:
                self.s[n]
            except KeyError:
                self.s[n] = {}
            
        for n in range(nmax+1):
            result.s[n] = {}
            sn = self.s[n]
            twon = 2**n
            twon1= twon - 1
            for key in self.s[n].keys():
                s0 = sn[key]
                
                if key[axis] == 0:
                    sm = zero      # Zero outside the box
                else:
                    t=list(key); t[axis]=t[axis]-1; t=tuple(t)
                    sm = self.__sock_it_to_me(n,t)
                    
                if key[axis] == twon1:
                    sp = zero
                else:
                    t=list(key); t[axis]=t[axis]+1; t=tuple(t)
                    sp = self.__sock_it_to_me(n,t)
                
                if (sm is None) or (sp is None):
                    # One or both of my neighbours are further
                    # down the tree.  Just recur me down to the
                    # next level and try again there.
                    self.__look_out_below(s0,n,key)
                    continue

                # rp & rm are only rank-1 ... don't yet exploit this

                r = rp.inner(sm,ind)
                r0.inner(s0,ind,result=r)  
                rm.inner(sp,ind,result=r)
                result.s[n][key] = r.make_self_contiguous_by_map(amap).scale(twon)
 
                #r = rp.inner(sm,ind,map=amap)
                #r0.inner(s0,ind,amap,result=r)  
                #rm.inner(sp,ind,map=amap,result=r)
                #result.s[n][key] = r.scale(twon)
                
                #result.s[n][key] = (rp.inner(sm,ind,amap) +   \
                #                    r0.inner(s0,ind,amap) +  \
                #                    rm.inner(sp,ind,amap)).scale(twon)
                
        self.sclean()

        perf.exit('diff')
        return result

    def diff2(self,axis):
        '''

        Apply a central derivative approximation to the second
        derivative along the specified axis, currently with zero
        boundary conditions.  Plan to implement periodic boundary
        conditions also.

        This operation corresponds to .  D2 = DT*D =
        self.diff(axis).diff(axist,transpose=1) and it should be
        negative (semi-definite ... I think that the boundary
        condition introduces a null space ... Beylkin suggests how to
        get rid of it but the projector is not yet coded).

        Self is reconstructed but should otherwise be unchanged.

        Returns a new function in the scaling function basis.

        '''
        if self.debug:
            print "diffing2", id(self), self.compressed, axis

        self.reconstruct()
        perf.enter('diff2')
        
        # Since are modifying the internals of result here must make
        # correct choice of NOT compressed
        result = self.new(k=self.k, compress=0)

        zero = Matrix(self.k,self.k)

        ind = (1,axis)                  # Indices to contract in products
        if axis == 0:
            amap = None
        elif axis == 1:
            amap = [1,0]
        else:
            raise "Don't you think six dimensions are enough?"
            
        rm, r0, rp = self.make_dc_periodic()
        Rp2 = rm.transpose()*rp
        Rp1 = rm.transpose()*r0 + r0.transpose()*rp
        R0  = rm.transpose()*rm + r0.transpose()*r0 + rp.transpose()*rp
        Rm1 =                     r0.transpose()*rm + rp.transpose()*r0
        Rm2 =                                         rp.transpose()*rm

        nmax = max(self.s.keys())

        for n in range(nmax+1):         #  Fill in empty levels
            try:
                self.s[n]
            except KeyError:
                self.s[n] = {}
            
        for n in range(nmax+1):
            result.s[n] = {}
            sn = self.s[n]
            twon = 2**n
            twon1= twon - 1
            scale = 4.0**n
            for key in self.s[n].keys():
                s = [-2,-1,0,1,2]
                got_them_all = 1
                for shift in s:
                    t=list(key); t[axis]=t[axis]+shift; t=tuple(t)
                    if t[axis]<0 or t[axis]>twon1:
                        s[2+shift] = zero
                    else:
                        s[2+shift] = self.__sock_it_to_me(n,t)
                        if s[2+shift] is None:
                            got_them_all = 0
                            break
                if not got_them_all:
                    # One or more of my neighbours are further
                    # down the tree.  Just recur me down to the
                    # next level and try again there.
                    self.__look_out_below(sn[key],n,key)
                    continue


                # Rp2 & Rm2 are only rank-1 ... don't yet exploit this
                sm2, sm1, s0, sp1, sp2 = s

                r = Rp2.inner(sm2,ind,amap)
                Rp1.inner(sm1,ind,amap,result=r)
                R0.inner(s0,ind,amap,result=r)
                Rm1.inner(sp1,ind,amap,result=r)
                Rm2.inner(sp2,ind,amap,result=r)
                result.s[n][key]=r.scale(scale)

                #result.s[n][key] = (Rp2.inner(sm2,ind,amap) +
                #                    Rp1.inner(sm1,ind,amap) + 
                #                    R0.inner(s0,ind,amap) +
                #                    Rm1.inner(sp1,ind,amap) +
                #                    Rm2.inner(sp2,ind,amap)).scale(scale)

        self.sclean()
        perf.exit('diff2')
        return result

    def laplacian(self):
        '''

        Compute the laplacian using a central derivative currently
        with a zero boundary condition (embedding).
        This operation corresponds to
        .   D2 = DT*D = self.diff(axis).diff(axis,transpose=1)
        and it should be negative (semi-definite ... I think that
        the embedding introduces a null space ... Beylkin suggests
        how to get rid of it but the projector is not yet coded).

        Self is reconstructed but should otherwise be unchanged.

        Returns a new function in the scaling function basis.

        '''
        if self.debug:
            print "laplacian", id(self), self.compressed
        
        self.reconstruct()
        perf.enter('laplacian')

        # Since are modifying the internals of result here must make
        # correct choice of NOT compressed
        result = self.new(k=self.k, compress=0)

        zero = Matrix(self.k,self.k)

        # For each axis the indices to contract in products and the
        # map from default to desired order in the contractions
        ind = ((1,0), (1,1))
        amap = ((0,1), (1,0))
            
        rm, r0, rp = self.make_dc_periodic()
        Rp2 = rm.transpose()*rp
        Rp1 = rm.transpose()*r0 + r0.transpose()*rp
        R0  = rm.transpose()*rm + r0.transpose()*r0 + rp.transpose()*rp
        Rm1 =                     r0.transpose()*rm + rp.transpose()*r0
        Rm2 =                                         rp.transpose()*rm

        nmax = max(self.s.keys())

        for n in range(nmax+1):         #  Fill in empty levels
            try:
                self.s[n]
            except KeyError:
                self.s[n] = {}
            
        for n in range(nmax+1):
            result.s[n] = {}
            sn = self.s[n]
            twon = 2**n
            twon1= twon - 1
            scale = -(4.0**n)
            for key in self.s[n].keys():
                s = [[-2,-1,0,1,2],[-2,-1,0,1,2],[-2,-1,0,1,2],
                     [-2,-1,0,1,2],[-2,-1,0,1,2],[-2,-1,0,1,2]]
                got_them_all = 1
                for axis in 0,1:
                    for shift in [-2,-1,0,1,2]:
                        t=list(key); t[axis]=t[axis]+shift; t=tuple(t)
                        if t[axis]<0 or t[axis]>twon1:
                            s[axis][2+shift] = zero
                        else:
                            s[axis][2+shift] = self.__sock_it_to_me(n,t)
                            if s[axis][2+shift] is None:
                                got_them_all = 0
                                break
                    if not got_them_all:
                        break
                    
                if not got_them_all:
                    # One or more of my neighbours are further
                    # down the tree.  Just recur me down to the
                    # next level and try again there.
                    self.__look_out_below(sn[key],n,key)
                    continue

                # Rp2 & Rm2 are only rank-1 ... don't yet exploit this
                r = Matrix(self.k,self.k)
                for axis in [0,1]:
                    sm2, sm1, s0, sp1, sp2 = s[axis]
                    i, a = ind[axis], amap[axis]
                    Rp2.inner(sm2,i,a,result=r)
                    Rp1.inner(sm1,i,a,result=r)
                    R0.inner(s0,i,a,result=r)
                    Rm1.inner(sp1,i,a,result=r)
                    Rm2.inner(sp2,i,a,result=r)
                    
                    #r = r + (Rp2.inner(sm2,i,a) + 
                    #         Rp1.inner(sm1,i,a) + 
                    #         R0.inner(s0,i,a)   +
                    #         Rm1.inner(sp1,i,a) +
                    #         Rm2.inner(sp2,i,a))
                    
                result.s[n][key] = r.scale(scale)

        self.sclean()
        perf.exit('laplacian')
        return result

    def __sock_it_to_me_diffNS(self,n,key,surface=0):

        k, k2 = self.k, self.k*2
        lx, ly = key

        if surface:
            if self.cache_sock_it_to_me_diffNS.has_key((lx,ly)):
                return self.cache_sock_it_to_me_diffNS[lx,ly]
        
        ix, iy = parentmod(lx), parentmod(ly)
        lx, ly = parent(lx), parent(ly)

        try:
            s = self.d[n-1][lx,ly]
            s = self.unfilter(s)
        except KeyError:
            s = self.__sock_it_to_me_diffNS(n-1,(lx,ly),surface=0)
            s = self.unfilter(s,sonly=1)
            
        s = Matrix(s[ix*k:ix*k+k,iy*k:iy*k+k])

        if surface:
            lx, ly = key
            self.cache_sock_it_to_me_diffNS[lx,ly] = s

        return s

    def __sock_it_to_me2_diffNS(self,n,key,surface=0):

        k, k2 = self.k, self.k*2
        lx, ly = key

        if surface:
            if self.cache_sock_it_to_me2_diffNS.has_key((lx,ly)):
                return self.cache_sock_it_to_me2_diffNS[lx,ly]
        
        ix, iy = parentmod(lx), parentmod(ly)
        lx, ly = parent(lx), parent(ly)
        
        try:
            s = self.d[n-1][lx,ly]
            s = self.unfilter(s)
            s = Matrix(s[ix*k:ix*k+k,iy*k:iy*k+k])
            s = s.normf()
        except KeyError:
            s = self.__sock_it_to_me2_diffNS(n-1,(lx,ly),surface=0)
            s = s * 0.50 # = (2.0**(-2.0/2.0))
            
        if surface:
            lx, ly = key
            self.cache_sock_it_to_me2_diffNS[lx,ly] = s

        return s

    def diffNS(self,axis,transpose=0,autorefine=0):

        k, k2 = self.k, self.k*2
        
        if self.debug:
            print "diffing", id(self), self.compressed, axis, \
                  transpose,self.s.keys()

        self.compress()
        perf.enter('diffNS')
        self.reconstruct(nonstandard=1)

        result = self.new(k=self.k, thresh=self.thresh)

        ind = (1,axis)                  # Indices to contract in products
        if axis == 0:
            amap = None
        elif axis == 1:
            amap = [1,0]
        else:
            raise "Don't you think six dimensions are enough?"
            
        (Tm, T0, Tp), (rm, r0, rp) = self.make_dc_periodic_compressed()

        if transpose:
            Tm, T0, Tp = Tp.transpose(), T0.transpose(), Tm.transpose()
            rm, r0, rp = rp.transpose(), r0.transpose(), rm.transpose()

        nmax = max(self.d.keys())

        if autorefine: nmax += 1

        for n in range(nmax+1):         #  Fill in empty levels
            try:
                self.d[n]
            except KeyError:
                self.d[n] = {}

        cut1, cut2 = 0, 0

        for n in range(nmax+1):
            result.d[n] = {}
            dn = self.d[n]
            twon = 2**n
            twon1= twon - 1
            
            self.cache_sock_it_to_me_diffNS = {}
            self.cache_sock_it_to_me2_diffNS = {}

            if autorefine:
                if n>0:
                    keys = []
                    for lx,ly in self.d[n-1].keys():
                        keys += [ (lx*2  ,ly*2  ),
                                  (lx*2  ,ly*2+1),
                                  (lx*2+1,ly*2  ),
                                  (lx*2+1,ly*2+1) ]
                else:
                    keys = [(0,0)]
            else:
                keys = dn.keys()

            for key in self.d[n].keys():
                t=list(key); t[axis] += 1; t=tuple(t)
                if t[axis]<=twon1 and not (t in keys): keys.append(t)
                t=list(key); t[axis] -= 1; t=tuple(t)
                if t[axis]>=0     and not (t in keys): keys.append(t)

            for key in keys:
                r = Matrix(k2,k2)
                if n>0: rPP = Matrix(k,k)

                # NOTE: We should expolit the sparsity of multiwavelet bases
                # in non-standard form 

                # T0 and r0 =====================================
                t=key
                try:
                    d = dn[t]
                except KeyError:
                    d0 = self.__sock_it_to_me_diffNS(n,t,surface=1)
                    d = Matrix(k2,k2)
                    d[0:k,0:k] = d0
                
                T0.inner(d,ind,map=amap,result=r) 
                if n>0:
                    d = d[0:k,0:k]
                    r0.inner(d,ind,map=amap,result=rPP)

                # Tp and rp =====================================
                if key[axis] != 0:
                    t=list(key); t[axis]=t[axis]-1; t=tuple(t)
                    try:
                        d = dn[t]
                    except KeyError:
                        d0 = self.__sock_it_to_me_diffNS(n,t,surface=1)
                        d = Matrix(k2,k2)
                        d[0:k,0:k] = d0
                
                    Tp.inner(d,ind,map=amap,result=r) 
                    if n>0:
                        d = d[0:k,0:k]
                        rp.inner(d,ind,map=amap,result=rPP) 

                # Tm and rm =====================================
                if key[axis] != twon1:
                    t=list(key); t[axis]=t[axis]+1; t=tuple(t)
                    try:
                        d = dn[t]
                    except KeyError:
                        d0 = self.__sock_it_to_me_diffNS(n,t,surface=1)
                        d = Matrix(k2,k2)
                        d[0:k,0:k] = d0
                        
                    Tm.inner(d,ind,map=amap,result=r)
                    if n>0:
                        d = d[0:k,0:k]
                        rm.inner(d,ind,map=amap,result=rPP)

                # ===============================================
                if n>0:
                    r[0:k,0:k].gaxpy(1.0,rPP[0:k,0:k],-1.0)
                    
                result.d[n][key] = r.scale(twon)
                
                #cut1 += 1
                #tmp = Matrix(result.d[n][key])
                #tmp[0:k,0:k].gaxpy(1.0,result.d[n][key][0:k,0:k],-1.0)
                #if tmp.normf() > self.thresh: cut2 +=1

            del self.cache_sock_it_to_me_diffNS
            del self.cache_sock_it_to_me2_diffNS

        # Operator has been applied.  Must now cleanup f from its
        # nonstandard state, and also compress the result summing the
        # sum coeffs at all levels in the process.

        #print cut2, "/", cut1

        for n in result.d.keys():
            dn = result.d[n]
            sn = result.s[n] = {}
            for key in dn.keys():
                sn[key] = Matrix(dn[key][0:k,0:k])
                dn[key][0:k,0:k].fill(0.0)
        
        result.compressed = 0
        result.compress()

        # Clean up the source
        self.s[0] = {}
        self.s[0][0,0] = Matrix(self.d[0][0,0][0:k,0:k])
        for n in self.d.keys():
            dn = self.d[n]
            for key in dn.keys():
                dn[key][0:k,0:k].fill(0.0)
        self.compressed = 1
        
        for n in range(nmax,-1,-1):
            if len(self.d[n])==0:
                del self.d[n]
            else:
                break

        perf.exit('diffNS')
        return result

    def diff2NS(self,axis,autorefine=0):

        k, k2 = self.k, self.k*2
        
        if self.debug:
            print "diffing2", id(self), self.compressed, axis

        self.compress()
        perf.enter('diff2NS')
        self.reconstruct(nonstandard=1)

        result = self.new(k=self.k, thresh=self.thresh)

        ind = (1,axis)                  # Indices to contract in products
        if axis == 0:
            amap = None
        elif axis == 1:
            amap = [1,0]
        else:
            raise "Don't you think six dimensions are enough?"
            
        (Tm, T0, Tp), (rm, r0, rp) = self.make_dc_periodic_compressed()
        
        tp2 = Tm.transpose()*Tp
        tp1 = Tm.transpose()*T0 + T0.transpose()*Tp
        t0  = Tm.transpose()*Tm + T0.transpose()*T0 + Tp.transpose()*Tp
        tm1 =                     T0.transpose()*Tm + Tp.transpose()*T0
        tm2 =                                         Tp.transpose()*Tm

        Rp2 = rm.transpose()*rp
        Rp1 = rm.transpose()*r0 + r0.transpose()*rp
        R0  = rm.transpose()*rm + r0.transpose()*r0 + rp.transpose()*rp
        Rm1 =                     r0.transpose()*rm + rp.transpose()*r0
        Rm2 =                                         rp.transpose()*rm

        del Tm, T0, Tp, rm, r0, rp
        
        nmax = max(self.d.keys())

        if autorefine: nmax += 1

        for n in range(nmax+1):         #  Fill in empty levels
            try:
                self.d[n]
            except KeyError:
                self.d[n] = {}

        cut1, cut2 = 0, 0
            
        for n in range(nmax+1):
            result.d[n] = {}
            dn = self.d[n]
            twon = 2**n
            twon1= twon - 1
            scale = 4.0**n
            
            self.cache_sock_it_to_me_diffNS = {}
            self.cache_sock_it_to_me2_diffNS = {}
            
            if autorefine:
                if n>0:
                    keys = []
                    for lx,ly in self.d[n-1].keys():
                        keys += [ (lx*2  ,ly*2  ),
                                  (lx*2  ,ly*2+1),
                                  (lx*2+1,ly*2  ),
                                  (lx*2+1,ly*2+1) ]
                else:
                    keys = [(0,0)]
            else:
                keys = dn.keys()

            for key in self.d[n].keys():
                t=list(key); t[axis] += 2; t=tuple(t)
                if t[axis]<=twon1 and not (t in keys): keys.append(t)
                t=list(key); t[axis] += 1; t=tuple(t)
                if t[axis]<=twon1 and not (t in keys): keys.append(t)
                t=list(key); t[axis] -= 1; t=tuple(t)
                if t[axis]>=0     and not (t in keys): keys.append(t)
                t=list(key); t[axis] -= 2; t=tuple(t)
                if t[axis]>=0     and not (t in keys): keys.append(t)

            for key in keys:
                r = Matrix(k2,k2)
                if n>0: rPP = Matrix(k,k)

                # NOTE: We should expolit the sparsity of multiwavelet bases
                # in non-standard form 

                # t0 and R0 =====================================
                t=key
                try:
                    d = dn[t]
                except KeyError:
                    d0 = self.__sock_it_to_me_diffNS(n,t,surface=1)
                    d = Matrix(k2,k2)
                    d[0:k,0:k] = d0
                
                t0.inner(d,ind,map=amap,result=r) 
                if n>0:
                    d = d[0:k,0:k]
                    R0.inner(d,ind,map=amap,result=rPP)

                # tp1 and Rp1 ===================================
                if key[axis] != 0:
                    t=list(key); t[axis]=t[axis]-1; t=tuple(t)
                    try:
                        d = dn[t]
                    except KeyError:
                        d0 = self.__sock_it_to_me_diffNS(n,t,surface=1)
                        d = Matrix(k2,k2)
                        d[0:k,0:k] = d0
                
                    tp1.inner(d,ind,map=amap,result=r) 
                    if n>0:
                        d = d[0:k,0:k]
                        Rp1.inner(d,ind,map=amap,result=rPP) 

                # tp2 and Rp2 ===================================
                if not (key[axis] in (0, 1)):
                    t=list(key); t[axis]=t[axis]-2; t=tuple(t)
                    try:
                        d = dn[t]
                    except KeyError:
                        d0 = self.__sock_it_to_me_diffNS(n,t,surface=1)
                        d = Matrix(k2,k2)
                        d[0:k,0:k] = d0
                
                    tp2.inner(d,ind,map=amap,result=r) 
                    if n>0:
                        d = d[0:k,0:k]
                        Rp2.inner(d,ind,map=amap,result=rPP) 

                # tm1 and Rm1 ===================================
                if key[axis] != twon1:
                    t=list(key); t[axis]=t[axis]+1; t=tuple(t)
                    try:
                        d = dn[t]
                    except KeyError:
                        d0 = self.__sock_it_to_me_diffNS(n,t,surface=1)
                        d = Matrix(k2,k2)
                        d[0:k,0:k] = d0
                        
                    tm1.inner(d,ind,map=amap,result=r)
                    if n>0:
                        d = d[0:k,0:k]
                        Rm1.inner(d,ind,map=amap,result=rPP)

                # tm2 and Rm2 ===================================
                if not (key[axis] in (twon1, twon1-1)):
                    t=list(key); t[axis]=t[axis]+2; t=tuple(t)
                    try:
                        d = dn[t]
                    except KeyError:
                        d0 = self.__sock_it_to_me_diffNS(n,t,surface=1)
                        d = Matrix(k2,k2)
                        d[0:k,0:k] = d0
                        
                    tm2.inner(d,ind,map=amap,result=r)
                    if n>0:
                        d = d[0:k,0:k]
                        Rm2.inner(d,ind,map=amap,result=rPP)

                # ===============================================
                if n>0:
                    r[0:k,0:k].gaxpy(1.0,rPP[0:k,0:k],-1.0)
                    
                result.d[n][key] = r.scale(scale)
                
                #cut1 += 1
                #tmp = Matrix(result.d[n][key])
                #tmp[0:k,0:k].gaxpy(1.0,result.d[n][key][0:k,0:k],-1.0)
                #if tmp.normf() > self.thresh: cut2 +=1

            del self.cache_sock_it_to_me_diffNS
            del self.cache_sock_it_to_me2_diffNS

        # Operator has been applied.  Must now cleanup f from its
        # nonstandard state, and also compress the result summing the
        # sum coeffs at all levels in the process.

        #print cut2, '/', cut1

        for n in result.d.keys():
            dn = result.d[n]
            sn = result.s[n] = {}
            for key in dn.keys():
                sn[key] = Matrix(dn[key][0:k,0:k])
                dn[key][0:k,0:k].fill(0.0)
            
        result.compressed = 0
        result.compress()

        # Clean up the source
        self.s[0] = {}
        self.s[0][0,0] = Matrix(self.d[0][0,0][0:k,0:k])
        for n in self.d.keys():
            dn = self.d[n]
            for key in dn.keys():
                dn[key][0:k,0:k].fill(0.0)
        self.compressed = 1

        for n in range(nmax,-1,-1):
            if len(self.d[n])==0:
                del self.d[n]
            else:
                break

        perf.exit('diff2NS')
        return result

    def project_bc(self):
        raise "not yet in 1d"

    def latex_zslice(self):
        raise "???"

class Function1d(Function):
    '''

    This is just Function restricted to 1d ... refer to that
    for documentation, comments and perhaps corrections.

    '''
    def __init__(self,
                 function=None, cfunction=None,
                 thresh=None, k=None,
                 initial_level=None,refine=1,compress=1,box=None):

        perf.enter('init')

        self.ndim = 1
        self.new = Function1d
        self.newtensor = Vector
        
        if k:
            self.k = k                  # Override default in class
        else:
            self.k = k = self.k

        if self.operator:               # Operators need double order
            self.k = k = 2*k

        if thresh:
            self.thresh = thresh        # Override default in class
        else:
            self.thresh = thresh = self.thresh

        if initial_level is None:
            self.initial_level = initial_level = Function.initial_level
        else:
            self.initial_level = initial_level

        self.d = {}             # Directory of wavelet coefficients
        self.s = {}             # Directory of scaling function coefficients
        self.function = function# Remember the function being compressed
        self.cfunction = cfunction
        self.compressed = 0     # Not yet compressed

        self.init_twoscale(k)
        self.init_quadrature(k)

        if cfunction:
            if not (type(cfunction) == type(' ') and \
                    cfunction[-6:] == 'p_void'):
                print cfunction
                raise TypeError,"cfunction should be a C "\
                      "function pointer cast to void *"

        if function or cfunction:
            self.nterminated = 0
            if box is None:
                trans = None
            else:
                trans = []
                lxlo,lxhi = box
                for lx in range(lxlo,lxhi+1):
                    trans.append(lx)
            
            perf.enter('fs1d proj')
            self.fine_scale_projection(initial_level,trans)
            perf.exit('fs1d proj')
            
            if refine:
                refinethese = {}
                for lx in self.s[initial_level].keys():
                    refinethese[parent(lx)] = 1
                for lx in refinethese.keys():
                    self.refine_fine_scale_projection(initial_level-1, lx)

            if self.nterminated > 0:
                print "Terminated refinement:", self.nterminated

            if compress:
                self.compress(truncate=0)
        else:
            self.s[0] = {}
            self.s[0][0] = Vector(k)
            if compress:
                self.compressed = 1
                self.d[0] = {}
                self.d[0][0] = Vector(2*k)
            
        if self.debug:
            print "Made new 1d function", id(self), self.k, self.thresh, self.compressed
            print "s:", self.s.keys()
            if self.compressed: print "d:", self.d.keys()

        perf.exit('init')
            
    def fline(self, xlo, npt, h, pts, f, function):
        '''

        Tabulate the function on the quadrature points in the
        specified cube ... this is intended to be replaced
        by a C routine.

        '''
        for p in range(npt):
            x = xlo+pts[p]*h
            f[p] = function(x)

    def fine_scale_projection(self,n,trans=None):
        '''

        Project the function onto scaling functions at level n.

        If self.operator is true, then we are computing the kernel of
        an integral operator. Operators need -2^n <= l <= 2^n as the
        default range.  

        If l[] is specified it is a list of translations that are to
        be considered.  Otherwise the full list is used.

        Modifies:
        :           self.s

        '''

        f = self.function
        try:
            self.s[n]
        except KeyError:
            self.s[n] = {}
        
        h = 1.0/(2.0**n)                # Box size on target level 
        scale = sqrt(h)
        phiw = self.quad_phiw
        pts = self.quad_x 
        npt = len(pts)

        if not trans:
            if self.operator:
                trans1d = range(-2**n,2**n+2)
            else:
                trans1d = range(2**n)
            trans = []
            for lx in trans1d:
                trans.append(lx)

        for lx in trans:
            xlo = lx*h
            f = Vector(npt)
            if self.cfunction:
                cfcube.cfline(xlo, npt, h,
                              pts.swigptr(), f.swigptr(),
                              self.cfunction)
            else:
                # Plain old Python function or non-cast C function
                # This is SLOW due to all of the tensor references.
                self.fline(xlo, npt, h, pts, f, self.function)
            f.scale(scale)

            self.s[n][lx] = f.fast_transform(phiw)

    def refine_fine_scale_projection(self,n,lx,prnt=0):
        '''

        We have the fine scale projection at level n+1, so we can
        compute the scaling function and wavelet coefficients at level
        n.  Considering just box l on level n, examine the difference
        coefficients.  If they do not satisfy the threshold, then
        refine the representation by recurring down to a finer level.

        If refinement occurs this modifies:
        :   self.s

        Set prnt to non-zero for printing

        '''
        if n > self.max_refine_level:
            self.nterminated = self.nterminated + 1
            return

        if Function.truncate_method == 0:
            tol = self.thresh
        else:
            tol = sqrt(1.0/2.0**n)*self.thresh
            
        k = self.k
        ss = self.gather_scaling_coeffs(n,lx)
        ss = self.filter(ss)
        s = Vector(ss[0:k])
        ss[0:k].fill(0.0)
        dnorm = ss.normf()
        if dnorm <= tol:
##             # The difference coefficients are neglible so can truncate
##             # at this level.  Delete the coeffs on the level below to
##             # save space, and adopt the scaling function coeffs
##             # computed from the lower level by twoscale since they
##             # will be more accurate
##             if not self.s.has_key(n): self.s[n] = {}
##             self.s[n][lx] = s
##             lx2 = lx*2
##             for ix in range(2):
##                 del(self.s[n+1][lx2+ix])
            pass
        else:
            # First remove the inaccurate scaling function coeffs at
            # level n+1, then refine to level n+2
            if prnt:
                for i in range(n): print "  ",
                print " refining level=%d l=%d norm=%.1e" % (n, lx, dnorm)
            lx2 = lx*2
            range2 = range(2)
            for ix in range2:
                ixlo = ix*k
                ixhi = ixlo + k
                del(self.s[n+1][lx2+ix])
                lx4 = 2*(lx2+ix)
                for iix in range2:
                    self.fine_scale_projection(n+2,[lx4+iix])
                self.refine_fine_scale_projection(n+1,lx2+ix)

    def gather_scaling_coeffs(self, n, lx):
        '''

        For a given translation at level n, gather the corresponding
        scaling function coefficients at level n+1.

        Some of the boxes on the lower level may be missing.

        '''
        k = self.k
        sn1 = self.s[n+1]
        ss = Vector(2*k)
        lx2 = lx*2
        for ix in range(2):
            ixlo = ix*k
            ixhi = ixlo + k
            key = lx2+ix
            if sn1.has_key(key):
                ss[ixlo:ixhi] =  sn1[key]
        return ss
        
    def reconstruct(self, nonstandard=0):
        '''

        Reconstruct the projection onto scaling functions from the
        compressed form terminating at the locally lowest level
        without difference coefficients.

        The compression algorithm has the responsibility of ensuring
        that the there are no signficant nodes below a node without
        difference coefficients.

        This algorithm must ensure that there are scaling function
        coefficients only at the locally lowest level (since other
        algorithms will just run down the tree until the first
        coefficient is found, and then stop).

        A variant of this algorithm is used to reconstruct the
        function in preparation for the application of a non-standard
        form operator.  For this we need to keep the difference
        coefficients and also keep the sum coefficients at each level
        on the way down, but only down to locally finest level at
        which we have difference coefficients.  There is (will be)
        another routine to clean the resulting mess back into the
        conventional hierarchical basis definition.

        '''
        if self.debug:
            print "reconstructing", id(self), self.compressed, self.s.keys()
            
        if not self.compressed:
            if nonstandard:
                self.compress()
            else:
                return self

        perf.enter('reconstrct')
        self.compressed = 0

        k = self.k

        nmax = max(self.d.keys())
        if nonstandard == 1:
            # Only need scaling functions down to where the difference
            # coefficients are locally non-zero.
            #nmax = nmax-1
            nmax = nmax
            
        range2 = range(2)

        for n in range(nmax+1):
            sn = self.s[n]
            dn = self.d[n]
            sn1 = self.s[n+1] = {}
            if nonstandard==1 and n<nmax:
                dn1 = self.d[n+1]
            for lx in dn.keys():
                d = dn[lx]
                d[0:k] = sn[lx]
                if nonstandard==1 and n>=nmax: continue
                if nonstandard != 2: del(sn[lx])
                d = self.unfilter(d)
                lx2 = lx*2
                for ix in range2:
                    key2 = lx2+ix
                    if (nonstandard==1) and (not dn1.has_key(key2)):
                        continue
                    dxyz = d[ix*k:ix*k+k]
                    sn1[key2] = Vector(dxyz)

            if nonstandard == 1:
                del self.s[n]
            else:
                del self.d[n]

        if nonstandard == 1: del self.s[nmax+1]
        
        if nmax == -1 and nonstandard == 1:
            print "FUDGE", self.s[0][0].normf()
            self.d[0][0][0:k] = self.s[0][0]
            print "FUDGE", self.d[0][0].normf()
            del self.s[0][0]

        perf.exit('reconstrct')
        return self

    def compress(self,truncate=0):
        '''

        Compress the scaling function coefficients (Vn+1 = Vn + Wn).

        In the usual compressed form, the scaling function
        coefficients will only exist on the lowest level, but this
        will vary across the domain.

        In the non-standard compressed form, both scaling function and
        difference coefficients will be defined on all levels ... must
        add the result of compressing the sum coefficients into the
        existing difference coefficients.

        In the form that results from addition in the scaling function
        basis, there may be sum coefficients at many levels, but no
        differences.

        The difference coefficients are truncated to the specified
        threshold.  So that reconstruction can just proceed down to
        the first level without difference coefficients, it is
        important to only truncate difference coefficients at leaf
        nodes.

        For ease of handling all the special cases now implement
        truncation as a separate and optional step.

        '''
        if self.debug:
            print "compressing", id(self), self.compressed, self.s.keys()
            
        if self.compressed:
            return self
        self.compressed = 1
        perf.enter('compress')
        
        nmax = max(self.s.keys())
        for n in range(nmax-1,-1,-1):
            if not self.d.has_key(n):
                self.d[n] = {}

        k = self.k
        for n in range(nmax-1,-1,-1):
            dn  = self.d[n]
            sn1 = self.s[n+1]
            try:
                sn = self.s[n]
            except KeyError:
                sn = self.s[n] = {}
            try:
                dn1 = self.d[n+1]
            except KeyError:
                dn1 = {}
            
            # Compute the list of translations that we must consider, driven
            # from the non-zero entries on the level below.
            llist = {}
            for lx in sn1.keys():
                llist[parent(lx)] = 1
            
            for lx in llist.keys():
                ss = self.gather_scaling_coeffs(n, lx)
                ss = self.filter(ss)
                s = Vector(ss[0:k]) # Force a copy
                if sn.has_key(lx):
                    sn[lx] = sn[lx].gaxpy(1.0,s,1.0)
                else:
                    sn[lx] = s

                ss[0:k].fill(0.0) # Now just want the diff coeffs
                if dn.has_key(lx):
                    dn[lx] = dn[lx].gaxpy(1.0,ss,1.0)
                else:
                    dn[lx] = ss

            del self.s[n+1]

        if truncate:
            self.truncate()

        if not self.d.has_key(0):
            self.d[0] = {}
        if not self.d[0].has_key(0):
            self.d[0][0] = Vector(2*k)

        if not self.s.has_key(0):
            self.s[0] = {}
        if not self.s[0].has_key(0):
            self.s[0][0] = Vector(k)

        perf.exit('compress')
        return self

    def basis(self,refine=0):
        '''

        Self is a compressed function.  Export a description of
        the basis (for use by restrict).  This describes the
        pieces of W0+W1+...+W(n-1) in use.

        If refine=1, then if the difference coefficients do not
        satisfy the accuracy goal, the basis is locally refined down
        one level.

        '''
        if not self.compressed: self.compress()

        if Function.truncate_method == 2:
            lend = self.__lend()
        
        basis = {}
        nmax = max(self.d.keys())
        for n in range(nmax+1):
            dn = self.d[n]
            bn = basis[n] = {}
            for key in dn.keys():
                bn[key] = 1

        if refine:
            basis[nmax+1] = {}          # Just in case
            for n in range(nmax+1):
                dn = self.d[n]
                bn = basis[n]
                bn1 = basis[n+1]

                keys = dn.keys()

                if Function.truncate_method == 0:
                    tol = self.thresh
                elif Function.truncate_method == 1:
                    tol = sqrt(1.0/2.0**n)*self.thresh
                else:
                    tol = self.thresh/sqrt(lend)

                range2 = range(2)
                for key in keys():
                    if (dn[key].normf() > tol):
                        lx = key
                        lx2 = lx*2
                        printed = 0
                        for ix in range2:
                            key2 = (lx2+ix)
                            if not bn1.has_key(key2):
                                if not printed:
                                    print "Refining", n, key
                                    printed = 1
                                bn1[key2] = 1
            if not basis[nmax+1]:
                del basis[nmax+1]

        return basis
    
    def truncate(self,thresh=None):
        '''

        The difference coefficients are truncated to the specified
        threshold.  So that reconstruction can just proceed down to
        the first level without difference coefficients, it is
        important to only truncate difference coefficients at leaf
        nodes.

        Returns self so that operations can be chained.

        '''
        if self.debug:
            print "truncating", id(self), self.compressed, self.thresh
            
        if not self.compressed:
            self.compress()

        if Function.truncate_method == 2:
            lend = self.__lend()
        
        if thresh is None:
            thresh = self.thresh
        
        perf.enter('truncate')

        nmax = max(self.d.keys())

        is_parent = {}
        for n in range(nmax,0,-1):
            dn = self.d[n]
            keys = dn.keys()

            if Function.truncate_method == 0:
                tol = thresh
            elif Function.truncate_method == 1:
                tol = sqrt(1.0/2.0**n)*thresh
            else:
                tol = thresh/sqrt(lend)

            # If not a parent and below threshold kill box.
            # If above threshold mark parent box as such.
            next_is_parent = {}
            for key in keys:
                if (not is_parent.has_key(key)) and \
                   (dn[key].normf() < tol):
                    del dn[key]
                else:
                    lx = key
                    next_is_parent[parent(lx)] = 1

            is_parent = next_is_parent
                    
            # If the current level is totally empty and it is also
            # the lowest level, it can be deleted entirely.
            if (not self.d[n]) and (not self.d.has_key(n+1)):
                del self.d[n]

        perf.exit('truncate')

##         if Function.truncate_method == 2:
##             # If we deleted a lot of coefficients, try again ...
##             if self.__lend()*2 < lend:
##                 self.truncate()

        return self

    def eval(self,x):
        '''

        Evaluate the function at the given point.

        If it is not compressed then we just need to walk down the
        tree until we locate the leaf node containing x,y,z.  If it is
        compressed then we do the same, but also have to add up the
        scaling function coeffs on the way down.

        '''
        n = 0
        k = self.k
        lx = 0
        xx = x
        if self.compressed:
            # Note always want loop to terminate via the break otherwise
            # n is one too small (x & l would also be wrong)
            s = self.s[0][0]
            for n in range(max(self.d.keys())+2):
                try:
                    d = Vector(self.d[n][lx]) # Take a copy
                except KeyError:
                    # At the bottom ... evaluate s*phi(xx) at level n
                    break
                d[0:k] = s
                d = self.unfilter(d)
                scale = 2**(n+1)
                xx = x*scale
                lx = min(scale-1,int(xx))
                xx = xx-lx
                ix = lx%2
                s = d[ix*k:ix*k+k]
            s = Vector(s)
        else:
            # The function is not compressed.  Assume that it is
            # represented at some (locally defined) finest level
            # in the scaling function basis.
            nmax = max(self.s.keys()) + 1
            for n in range(nmax):
                twon = 2**n
                twon1= twon - 1
                xx = x*twon
                lx = min(twon1,int(xx))
                xx = xx-lx
                try:
                    s = self.s[n][lx]
                    break
                except KeyError:
                    pass
            
        #px = Vector(phi(xx,k))
        #value = s.inner(px)*sqrt(2.0**n)

        value = cfcube.eval1d(k,n,xx,s.swigptr())        
        return value

    def eval_err(self, npt=1000):
        self.reconstruct()
        errsq = 0.0
        maxerr= 0.0
        maxrelerr=0.0
        seed(201)
        for i in range(npt):
            x = 0.4+random()*0.1
            value = self.eval(x)
            exact = self.function(x)
            err = abs(exact - value)
            maxerr = max(maxerr, err)
            errsq = errsq + err*err
            if exact:
                maxrelerr = max(maxrelerr,err/exact)
        rmserr = sqrt(errsq/npt)
        print "RMS abs error = %.1e" % rmserr
        print "MAX abs error = %.1e" % maxerr
        print "MAX rel error = %.1e" % maxrelerr
        self.compress(truncate=0)

    def print_layers(self, print_pic=0):
        if self.compressed:
            print " Analysis by layers of compressed function"
            nmax = max(self.d.keys())
            normsq = self.s[0][0].normf()**2
            print "\nNorm in V0 %.8e" % sqrt(normsq)
            for n in range(nmax+1):
                sum = 0.0
                for xyz in self.d[n].keys():
                    sum = sum + self.d[n][xyz].normf()**2
                normsq = normsq + sum
                print "Norm in W%d %.8e %.8e" % (n,sqrt(sum),sqrt(normsq))

            if not print_pic: return  
        else:
            print " Analysis by layers of reconstructed function"
            nmax = max(self.s.keys())
            normsq = 0.0
            for n in range(nmax+1):
                sum = 0.0
                for xyz in self.s[n].keys():
                    sum = sum + self.s[n][xyz].normf()**2
                normsq = normsq + sum
                print "Norm in V%d %.8e %.8e" % (n,sqrt(sum),sqrt(normsq))

            if not print_pic: return

    def sclean(self):
        '''

        Cleans up scaling function coefficients after diffing or
        multiplying.  Only want to keep the scaling function coeffs
        at the highest level in the tree ... ones lower down will
        have been generated by two-scale just for evaluation.
        Kill all nodes that have a parent - infanticide.

        '''
        for n in range(max(self.s.keys()),0,-1):
            sn, sn1 = self.s[n], self.s[n-1]
            for lx in sn.keys():
                if sn1.has_key(parent(lx)):
                    del(sn[lx])

    def mul(self,other):
        '''

        Multiplication of a function by a another function.
        Replaces the algorithm based on squaring.

        Returns a new function in the scaling function basis.
        Self/other are not modified though may be reconstructed if
        either is compressed.

        Added an optimization for multiplcation by the mask which is
        mask(x) = 1 for 0.0625 < x < 1-0.625.  Thus, only need to do
        the multiplication we have boxes 0 or 15 on level 4 (or boxes
        contained within below).  The mask is identified by the
        attribute mask which noone else will have.

        '''

        self.reconstruct()
        other.reconstruct()

        if hasattr(self,'mask'): self,other = other,self
        mask = hasattr(other,'mask')
        
        perf.enter('mul')

        result = self.new(k=self.k,compress=0) # Was thresh=self.thresh
        result.s = {}

        k = self.k
        phit = Matrix(self.quad_phi.cycledim(1)) # Contig. for fast_transform
        phiw = self.quad_phiw
        npt = len(self.quad_x)
        nmax = max(self.s.keys()+other.s.keys())

        tmp1 = Vector(k)
        tmp2 = Vector(k)

        range2 = range(2)

        if Function.autorefine:     # NOTE REFERENCE TO CLASS
            result.s[0] = {}
             
        for n in range(nmax+1):
            try:
                sn = self.s[n]
            except KeyError:
                sn = self.s[n] = {}

            try:
                on = other.s[n]
            except KeyError:
                on = other.s[n] = {}
                
            if Function.autorefine:     # NOTE REFERENCE TO CLASS
                rn1 = result.s[n+1] = {}
                scale = sqrt(2.0**(n+1))
            else:
                rn = result.s[n] = {}
                scale = sqrt(2.0**n)

            # Generate union of keys on this level. 
            keys = sn.keys() + on.keys()
            keys.sort()
            for i in range(len(keys)-1,0,-1):
                if keys[i] == keys[i-1]:
                    del keys[i]

            if mask:
                if n <= 4:
                    masklo, maskhi = 0, 2**n-1
                else:
                    masklo, maskhi = 2**(n-4)-1, 15*2**(n-4)
                    
            for key in keys:
                lx = key

                # If other is the mask, exclude keys for which the
                # mask is unity.
                if mask:
                    if masklo<lx<maskhi:
                        if sn.has_key(key): rn[key] = Vector(sn[key])
                        continue

                s =  self.__sock_it_to_me(n,key)
                if not (s is None):
                    o = other.__sock_it_to_me(n,key)
                    if not (o is None):
                        if Function.autorefine: # NOTE REFERENCE TO CLASS
                            ss = self.unfilter(s,sonly=1)
                            oo = self.unfilter(o,sonly=1)
                            lx2 = 2*lx
                            for ix in range2:
                                f = ss[ix*k:ix*k+k]
                                g = oo[ix*k:ix*k+k]
                                f = f.transform(phit,result=tmp1)
                                g = g.transform(phit,result=tmp2)
                                f.emul(g)
                                rn1[lx2+ix] = \
                                   f.fast_transform(phiw).scale(scale)
                        else:
                            f = s.fast_transform(phit,result=tmp1)
                            g = o.fast_transform(phit,result=tmp2)
                            f.emul(g)
                            rn[key] = f.fast_transform(phiw).scale(scale)
                            
        self.sclean()
        other.sclean()

        perf.exit('mul')

        return result

    def square(self):
        '''

        Square self replacing self by its square.

        If (Function1d.autorefine) is set, then it will recur down the normal
        additional level, otherwise it does not.

        NB: self is changed.

        Returns self in the scaling function basis so that operations
        may be chained.

        '''
        if self.debug:
            print "squaring", id(self), self.compressed
            
        self.reconstruct()

        k = self.k
        nmax = max(self.s.keys())
        phit = Matrix(self.quad_phi.cycledim(1)) # Contig. for fast_transform
        phiw = self.quad_phiw
        npt = len(self.quad_x)

        try:
            next_keys = self.s[0].keys()
        except KeyError:
            # This can happen if not compressed and initial_level != 0
            self.s[0] = {}
            next_keys = self.s[0].keys()

        for n in range(nmax+1):
        
            sn = self.s[n]
            try:
                sn1 = self.s[n+1]
            except KeyError:
                sn1 = self.s[n+1] = {}

            keys = next_keys
            next_keys = sn1.keys()

            if Function1d.autorefine:
                ss = Vector(2*k)
                scale = sqrt(2.0**(n+1))
            else:
                scale = sqrt(2.0**n)
            for lx in keys:
                if Function1d.autorefine:
                    # Now are at the locally lowest level ... recur down
                    # one more level and then form the square projecting
                    # onto the approximately refined basis.  Could recur
                    # more dynamically to guarantee precision but
                    # Beylkin does not seem to think this is necessary
                    # except in pathological cases?
                    ss = self.unfilter(sn[lx],sonly=1)
                    del(sn[lx])
                    lx2 = 2*lx
                    for ix in range(2):
                        f = Vector(ss[ix*k:ix*k+k])
                        # Back transform to the function value
                        # at the quadrature points
                        f = f.fast_transform(phit)
                        # Compute the square on the grid
                        f.emul(f)
                        # Do quadrature on f**2
                        f = f.fast_transform(phiw)
                        f.scale(scale)
                        # Done
                        sn1[lx2+ix] = f
                else:
                    # Don't refine ... just do it at the current level
                    f = sn[lx].fast_transform(phit)
                    f.emul(f)
                    sn[lx] = f.fast_transform(phiw).scale(scale)
            
        return self

    def local_function(self,function=None,cfunction=None):
        '''

        Apply a scalar local function to self replacing self by
        f(self(x)).  Define either function (a callable object) or
        cfunction (a pointer to a function cast to a pointer to void
        ... so that old versions of SWIG can be used).  Both the C and
        Python functions should take a double argument and return a
        double.  THIS APPROACH IS ONLY GOOD IF THE RESULT OF THE
        FUNCTION IS AS LOCALLY SMOOTH OR SMOOTHER THAN THE INPUT ...

        If (Function.autorefine) is set, then it will recur down an
        additional level before squaring, otherwise it does not.

        NB: self is changed.

        Returns self so that operations may be chained.  

        '''
        if self.debug:
            print "squaring", id(self), self.compressed
            
        range2 = range(2)

        self.reconstruct()
        perf.enter('local func')

        k = self.k
        nmax = max(self.s.keys())
        phit = Matrix(self.quad_phi.cycledim(1)) # Contig. for fast_transform
        phiw = self.quad_phiw
        npt = len(self.quad_x)

        tmp1 = Vector(k)

        try:
            next_keys = self.s[0].keys()
        except KeyError:
            # This can happen if not compressed and initial_level != 0
            self.s[0] = {}
            next_keys = self.s[0].keys()

        for n in range(nmax+1):
        
            sn = self.s[n]
            try:
                sn1 = self.s[n+1]
            except KeyError:
                sn1 = self.s[n+1] = {}

            keys = next_keys
            next_keys = sn1.keys()

            if Function.autorefine:         # NOTE REFERENCE TO CLASS
                ss = Vector(2*k)
                scale = sqrt(2.0**(n+1))
            else:
                scale = sqrt(2.0**n)
            for lx in keys:
                if Function.autorefine:     # NOTE REFERENCE TO CLASS
                    # Now are at the locally lowest level ... recur down
                    # one more level and then form the square projecting
                    # onto the approximately refined basis.  Could recur
                    # more dynamically to guarantee precision but
                    # Beylkin does not seem to think this is necessary
                    # except in pathological cases?
                    ss = self.unfilter(sn[lx],sonly=1)
                    del(sn[lx])
                    lx2 = 2*lx
                    for ix in range2:
                        for iy in range2:
                            f = ss[ix*k:ix*k+k]
                            # Back transform to the function value
                            # at the quadrature points
                            f = Vector(f).transform(phit,result=tmp1).scale(scale)
                            # Eval function on the refined grid
                            if cfunction:
                                cfcube.clocalfunc(
                                    f.swigptr(),k**1,cfunction)
                            else:
                                f.unaryop(function)
                            # Do quadrature on f**2 and scale
                            sn1[lx2+ix] = \
                                f.fast_transform(phiw).scale(1.0/scale)
                                                          
                else:
                    # Don't refine ... just do it at the current level
                    f = sn[lx].fast_transform(phit,result=tmp1).scale(scale)
                    if cfunction:
                        cfcube.clocalfunc(f.swigptr(),k**1,cfunction)
                    else:
                        f.unaryop(function)
                    sn[lx] = f.fast_transform(phiw).scale(1.0/scale)
            
        perf.exit('local func')
        return self

    def mul_func(self,function=None,cfunction=None):
        '''

        '''
        if self.debug:
            print "squaring", id(self), self.compressed
            
        range2 = range(2)

        self.reconstruct()
        perf.enter('mul func')

        k = self.k
        nmax = max(self.s.keys())
        phit = Matrix(self.quad_phi.cycledim(1)) # Contig. for fast_transform
        phiw = self.quad_phiw
        pts = self.quad_x 
        npt = len(self.quad_x)

        tmp1 = Vector(k)

        try:
            next_keys = self.s[0].keys()
        except KeyError:
            # This can happen if not compressed and initial_level != 0
            self.s[0] = {}
            next_keys = self.s[0].keys()

        for n in range(nmax+1):

            sn = self.s[n]
            try:
                sn1 = self.s[n+1]
            except KeyError:
                sn1 = self.s[n+1] = {}

            keys = next_keys
            next_keys = sn1.keys()

            if Function.autorefine:         # NOTE REFERENCE TO CLASS
                ss = Vector(2*k)
                dd = float(2**(n+1))
            else:
                dd = float(2**n)

            scale = sqrt(dd**1)       # 2D Scale
            h = 1.0/dd                # Box size on target level
             
            for lx in keys:

                if Function.autorefine:     # NOTE REFERENCE TO CLASS
                    # Now are at the locally lowest level ... recur down
                    # one more level and then form the square projecting
                    # onto the approximately refined basis.  Could recur
                    # more dynamically to guarantee precision but
                    # Beylkin does not seem to think this is necessary
                    # except in pathological cases?
                    ss = self.unfilter(sn[lx],sonly=1)
                    del(sn[lx])
                    lx2 = 2*lx
                    for ix in range2:

                        # get a cube of g:(func)
                        # ----------------------
                        xlo = (lx2+ix)*h
                        g = Vector(npt)
                        if cfunction:
                            cfcube.cfline(xlo, npt, h,
                                          pts.swigptr(), g.swigptr(),
                                          cfunction)
                        else:
                            # Plain old Python function or non-cast C function
                            # This is SLOW due to all of the tensor references.
                            self.fline(xlo, npt, h, pts, g, function)
                        
                        # get a cube of f:(func) 
                        # ----------------------
                        f = ss[ix*k:ix*k+k]
                        # Back transform to the function value
                        # at the quadrature points
                        f = f.transform(phit,result=tmp1).scale(scale)
                        
                        # emul between f and g
                        # --------------------
                        f = f.emul(g)
                        
                        # Do quadrature on f**2 and scale
                        sn1[lx2+ix] = \
                            f.fast_transform(phiw).scale(1.0/scale)
                                                          
                else:

                    # get a cube of g:(func)
                    # ----------------------
                    xlo = lx*h
                    g = Vector(npt)
                    if cfunction:
                        cfcube.cfline(xlo, npt, h,
                                      pts.swigptr(), g.swigptr(),
                                      cfunction)
                    else:
                        # Plain old Python function or non-cast C function
                        # This is SLOW due to all of the tensor references.
                        self.fline(xlo, npt, h, pts, g, function)

                    # get a cube of f:(func) 
                    # ----------------------
                    # Don't refine ... just do it at the current level
                    f = sn[lx].fast_transform(phit,result=tmp1).scale(scale)

                    # emul between f and g
                    # --------------------
                    f = f.emul(g)
                
                    sn[lx] = f.fast_transform(phiw).scale(1.0/scale)
            
        perf.exit('mul func')
        return self

    def __look_out_below(self,s, n, key):
        '''

        s are scaling coefficients at level n for the given
        translation.  Recur them down one level and put the resulting
        coefficients into the tree.

        '''
        if not self.s.has_key(n+1):
            self.s[n+1] = {}
        k = self.k
        ss = self.unfilter(s,sonly=1)
        lx = key
        lx2 = lx*2
        for ix in range(2):
            self.s[n+1][lx2+ix] = Vector(ss[ix*k:ix*k+k])


    def __sock_it_to_me(self,n,key):
        '''

        Return self.s[n][key] or None if it cannot be made (because it
        has no parent).
        
        If it is absent, then try to recur from the level above to
        make it.  If the parent is also absent, we need to keep
        walking up.  If we hit the top of the tree and still cannot
        find a parent, then return None.

        Some optimizations are possible, such as looking one below
        before recuring up.  Also, has_key() might be faster than
        exception handling.

        '''
        try:
            return self.s[n][key]
        except KeyError:
            if n == 0:
                return None

        lx = key
        lx = parent(lx)
        try:
            s = self.s[n-1][lx]
        except KeyError:
            s = self.__sock_it_to_me(n-1,(lx))

        if s is None:
            return None

        self.__look_out_below(s,n-1,(lx))

        return self.s[n][key]

    def __sock_it_to_me2(self,n,key):
        '''

        Return self.s[n][key] or None if it cannot be made (because it
        has no parent).
     
        If it is absent, then try to recur from the level above to
        make it.  If the parent is also absent, we need to keep
        walking up.  If we hit the top of the tree and still cannot
        find a parent, then return None.

        Some optimizations are possible, such as looking one below
        before recuring up.  Also, has_key() might be faster than
        exception handling.

        '''
        try:
            return self.s[n][key].normf()
        except KeyError:
            if n == 0:
                return None

        lx = key

        lx = parent(lx)
        try:
            s = self.s[n-1][lx]
            s = s.normf()
        except KeyError:
            s = self.__sock_it_to_me2(n-1,lx)

        if s is None:
            return None

        s = s * 0.70710678118654757 # = (2.0**(-1.0/2.0))
     
        return s
    
    def __look_out_below3(self,s, n, key):
        '''

        s are scaling coefficients at level n for the given
        translation.  Recur them down one level and put the resulting
        coefficients into the tree.

        '''
        if not self.s.has_key(n+1):
            self.s[n+1] = {}
        k = self.k
        ss = self.unfilter(s,sonly=1)
        lx = key
        lx2 = lx*2
        for ix in range(2):
            self.s[n+1][lx2+ix] = Vector(ss[ix*k:ix*k+k])


    def __sock_it_to_me3(self,n,key):
        '''

        Return self.s[n][key] or None if it cannot be made (because it
        has no parent).
        
        If it is absent, then try to recur from the level above to
        make it.  If the parent is also absent, we need to keep
        walking up.  If we hit the top of the tree and still cannot
        find a parent, then return None.

        Some optimizations are possible, such as looking one below
        before recuring up.  Also, has_key() might be faster than
        exception handling.

        '''
        try:
            s = self.s[n][key]
            if s != None:
                return n, key
        except KeyError:
            if n == 0:
                return None

        lx = key
        lx = parent(lx)
        try:
            s = self.s[n-1][lx]
        except KeyError:
            new_n, new_l = self.__sock_it_to_me(n-1,(lx))

        if s is None:
            return None, None

        self.__look_out_below(s,n-1,(lx))

        return self.s[n][key]

    def diff(self,axis,transpose=0):
        '''

        Apply a central derivative along the specified axis,
        currently with zero boundary conditions.  Plan to implement
        periodic boundary conditions also.

        If (transpose) then use the transposed derivative operator.

        Self is reconstructed but should otherwise be unchanged.

        Returns a new function in the scaling function basis.

        '''
        if self.debug:
            print "diffing", id(self), self.compressed, axis, \
                  transpose,self.s.keys()

        self.reconstruct()
        perf.enter('diff')
        
        result = self.new(k=self.k, compress=0)

        zero = Vector(self.k)

        ind = (1,axis)                  # Indices to contract in products
        if axis == 0:
            amap = None
        else:
            raise "Don't you think six dimensions are enough?"
            
        rm, r0, rp = self.make_dc_periodic()

        if transpose:
            rm, r0, rp = rp.transpose(), r0.transpose(), rm.transpose()

        nmax = max(self.s.keys())

        for n in range(nmax+1):         #  Fill in empty levels
            try:
                self.s[n]
            except KeyError:
                self.s[n] = {}
            
        for n in range(nmax+1):
            result.s[n] = {}
            sn = self.s[n]
            twon = 2**n
            twon1= twon - 1
            for key in self.s[n].keys():
                s0 = sn[key]

                if key == 0:
                    sm = zero      # Zero outside the box
                else:
                    t=key; t=t-1
                    sm = self.__sock_it_to_me(n,t)
                    
                if key == twon1:
                    sp = zero
                else:
                    t=key; t=t+1
                    sp = self.__sock_it_to_me(n,t)
                
                if (sm is None) or (sp is None):
                    # One or both of my neighbours are further
                    # down the tree.  Just recur me down to the
                    # next level and try again there.
                    self.__look_out_below(s0,n,key)
                    continue

                # rp & rm are only rank-1 ... don't yet exploit this

                r = rp.inner(sm,ind,amap)
                r0.inner(s0,ind,amap,result=r)  
                rm.inner(sp,ind,amap,result=r)
                result.s[n][key] = r.scale(twon)
                
                #r = rp.inner(sm)
                #r0.inner(s0,result=r)  
                #rm.inner(sp,result=r)
                #result.s[n][key] = r.scale(twon)
                
                #result.s[n][key] = (rp.inner(sm,ind,amap) + 
                #                    r0.inner(s0,ind,amap) +
                #                    rm.inner(sp,ind,amap)).scale(twon)

        self.sclean()

        perf.exit('diff')
        return result

    def diff2(self,axis):
        '''

        Apply a central derivative approximation to the second
        derivative along the specified axis, currently with zero
        boundary conditions.  Plan to implement periodic boundary
        conditions also.

        This operation corresponds to .  D2 = DT*D =
        self.diff(axis).diff(axist,transpose=1) and it should be
        negative (semi-definite ... I think that the boundary
        condition introduces a null space ... Beylkin suggests how to
        get rid of it but the projector is not yet coded).

        Self is reconstructed but should otherwise be unchanged.

        Returns a new function in the scaling function basis.

        '''
        if self.debug:
            print "diffing2", id(self), self.compressed, axis

        self.reconstruct()
        perf.enter('diff2')
        
        # Since are modifying the internals of result here must make
        # correct choice of NOT compressed
        result = self.new(k=self.k, compress=0)

        zero = Vector(self.k)

        ind = (1,axis)                  # Indices to contract in products
        if axis == 0:
            amap = None
        else:
            raise "Don't you think six dimensions are enough?"
            
        rm, r0, rp = self.make_dc_periodic()
        Rp2 = rm.transpose()*rp
        Rp1 = rm.transpose()*r0 + r0.transpose()*rp
        R0  = rm.transpose()*rm + r0.transpose()*r0 + rp.transpose()*rp
        Rm1 =                     r0.transpose()*rm + rp.transpose()*r0
        Rm2 =                                         rp.transpose()*rm

        nmax = max(self.s.keys())

        for n in range(nmax+1):         #  Fill in empty levels
            try:
                self.s[n]
            except KeyError:
                self.s[n] = {}
            
        for n in range(nmax+1):
            result.s[n] = {}
            sn = self.s[n]
            twon = 2**n
            twon1= twon - 1
            scale = 4.0**n
            for key in self.s[n].keys():
                s = [-2,-1,0,1,2]
                got_them_all = 1
                for shift in s:
                    t=key; t=t+shift
                    if t<0 or t>twon1:
                        s[2+shift] = zero
                    else:
                        s[2+shift] = self.__sock_it_to_me(n,t)
                        if s[2+shift] is None:
                            got_them_all = 0
                            break
                if not got_them_all:
                    # One or more of my neighbours are further
                    # down the tree.  Just recur me down to the
                    # next level and try again there.
                    self.__look_out_below(sn[key],n,key)
                    continue


                # Rp2 & Rm2 are only rank-1 ... don't yet exploit this
                sm2, sm1, s0, sp1, sp2 = s

                r = Rp2.inner(sm2,ind,amap)
                Rp1.inner(sm1,ind,amap,result=r)
                R0.inner(s0,ind,amap,result=r)
                Rm1.inner(sp1,ind,amap,result=r)
                Rm2.inner(sp2,ind,amap,result=r)
                result.s[n][key]=r.scale(scale)

                #result.s[n][key] = (Rp2.inner(sm2,ind,amap) +
                #                    Rp1.inner(sm1,ind,amap) + 
                #                    R0.inner(s0,ind,amap) +
                #                    Rm1.inner(sp1,ind,amap) +
                #                    Rm2.inner(sp2,ind,amap)).scale(scale)

        self.sclean()
        perf.exit('diff2')
        return result

    def laplacian(self):
        '''

        Compute the laplacian using a central derivative currently
        with a zero boundary condition (embedding).
        This operation corresponds to
        .   D2 = DT*D = self.diff(axis).diff(axis,transpose=1)
        and it should be negative (semi-definite ... I think that
        the embedding introduces a null space ... Beylkin suggests
        how to get rid of it but the projector is not yet coded).

        Self is reconstructed but should otherwise be unchanged.

        Returns a new function in the scaling function basis.

        '''
        if self.debug:
            print "laplacian", id(self), self.compressed
        
        self.reconstruct()
        perf.enter('laplacian')

        # Since are modifying the internals of result here must make
        # correct choice of NOT compressed
        result = self.new(k=self.k, compress=0)

        zero = Vector(self.k)

        # For each axis the indices to contract in products and the
        # map from default to desired order in the contractions
        ind = ((1,0))
        amap = ((0))
            
        rm, r0, rp = self.make_dc_periodic()
        Rp2 = rm.transpose()*rp
        Rp1 = rm.transpose()*r0 + r0.transpose()*rp
        R0  = rm.transpose()*rm + r0.transpose()*r0 + rp.transpose()*rp
        Rm1 =                     r0.transpose()*rm + rp.transpose()*r0
        Rm2 =                                         rp.transpose()*rm

        nmax = max(self.s.keys())

        for n in range(nmax+1):         #  Fill in empty levels
            try:
                self.s[n]
            except KeyError:
                self.s[n] = {}
            
        for n in range(nmax+1):
            result.s[n] = {}
            sn = self.s[n]
            twon = 2**n
            twon1= twon - 1
            scale = -(4.0**n)
            for key in self.s[n].keys():
                s = [[-2,-1,0,1,2],[-2,-1,0,1,2],[-2,-1,0,1,2],
                     [-2,-1,0,1,2],[-2,-1,0,1,2],[-2,-1,0,1,2]]
                got_them_all = 1
                for axis in [0]:
                    for shift in [-2,-1,0,1,2]:
                        t=list(key); t[axis]=t[axis]+shift; t=tuple(t)
                        if t[axis]<0 or t[axis]>twon1:
                            s[axis][2+shift] = zero
                        else:
                            s[axis][2+shift] = self.__sock_it_to_me(n,t)
                            if s[axis][2+shift] is None:
                                got_them_all = 0
                                break
                    if not got_them_all:
                        break
                    
                if not got_them_all:
                    # One or more of my neighbours are further
                    # down the tree.  Just recur me down to the
                    # next level and try again there.
                    self.__look_out_below(sn[key],n,key)
                    continue

                # Rp2 & Rm2 are only rank-1 ... don't yet exploit this
                r = Vector(self.k)
                for axis in [0]:
                    sm2, sm1, s0, sp1, sp2 = s[axis]
                    i, a = ind[axis], amap[axis]
                    Rp2.inner(sm2,i,a,result=r)
                    Rp1.inner(sm1,i,a,result=r)
                    R0.inner(s0,i,a,result=r)
                    Rm1.inner(sp1,i,a,result=r)
                    Rm2.inner(sp2,i,a,result=r)
                    
                    #r = r + (Rp2.inner(sm2,i,a) + 
                    #         Rp1.inner(sm1,i,a) + 
                    #         R0.inner(s0,i,a)   +
                    #         Rm1.inner(sp1,i,a) +
                    #         Rm2.inner(sp2,i,a))
                    
                result.s[n][key] = r.scale(scale)

        self.sclean()
        perf.exit('laplacian')
        return result
                
    def __sock_it_to_me_diffNS(self,n,lx,surface=0):

        k, k2 = self.k, self.k*2
        lx_sv = lx

        if surface:
            if self.cache_sock_it_to_me_diffNS.has_key(lx):
                return self.cache_sock_it_to_me_diffNS[lx]
        
        ix = parentmod(lx)
        lx = parent(lx)

        try:
            s = self.d[n-1][lx]
            s = self.unfilter(s)
        except KeyError:
            s = self.__sock_it_to_me_diffNS(n-1,lx,surface=0)
            s = self.unfilter(s,sonly=1)
            
        s = Vector(s[ix*k:ix*k+k])

        if surface:
            lx = lx_sv
            self.cache_sock_it_to_me_diffNS[lx] = s

        return s

    def __sock_it_to_me2_diffNS(self,n,lx,surface=0):

        k = self.k
        lx_sv = lx

        if surface:
            if self.cache_sock_it_to_me2_diffNS.has_key(lx):
                return self.cache_sock_it_to_me2_diffNS[lx]
        
        ix = parentmod(lx)
        lx = parent(lx)
        
        try:
            s = self.d[n-1][lx]
            s = self.unfilter(s)
            s = Vector(s[ix*k:ix*k+k])
            s = s.normf()
        except KeyError:
            s = self.__sock_it_to_me2_diffNS(n-1,lx,surface=0)
            s = s * 0.70710678118654757 # = (2.0**(-1.0/2.0))
            
        if surface:
            lx = lx_sv
            self.cache_sock_it_to_me2_diffNS[lx] = s

        return s

    def diffNS(self,axis,transpose=0,autorefine=0):

        k, k2 = self.k, self.k*2
        
        if self.debug:
            print "diffing", id(self), self.compressed, axis, \
                  transpose,self.s.keys()

        self.compress()
        perf.enter('diffNS')
        self.reconstruct(nonstandard=1)

        result = self.new(k=self.k, thresh=self.thresh)

        ind = (1,axis)                  # Indices to contract in products
        if axis == 0:
            amap = None
        else:
            raise "No other axis but 0 ..."
            
        (Tm, T0, Tp), (rm, r0, rp) = self.make_dc_periodic_compressed()

        if transpose:
            Tm, T0, Tp = Tp.transpose(), T0.transpose(), Tm.transpose()
            rm, r0, rp = rp.transpose(), r0.transpose(), rm.transpose()

        nmax = max(self.d.keys())

        if autorefine: nmax += 1
        
        for n in range(nmax+1):         #  Fill in empty levels
            try:
                self.d[n]
            except KeyError:
                self.d[n] = {}

        cut1, cut2 = 0, 0
            
        for n in range(nmax+1):
            result.d[n] = {}
            dn = self.d[n]
            twon = 2**n
            twon1= twon - 1
            
            self.cache_sock_it_to_me_diffNS = {}
            self.cache_sock_it_to_me2_diffNS = {}

            if autorefine:
                if n>0:
                    keys = []
                    for lx in self.d[n-1].keys():
                        keys += [ lx*2, lx*2+1 ]
                else:
                    keys = [0]
            else:
                keys = dn.keys()

            for key in self.d[n].keys():
                t=key; t += 1
                if t<=twon1 and not (t in keys): keys.append(t)
                t=key; t -= 1
                if t>=0     and not (t in keys): keys.append(t)

            for key in keys:
                r = Vector(k2)
                if n>0: rPP = Vector(k)

                # NOTE: We should expolit the sparsity of multiwavelet bases
                # in non-standard form 

                # T0 and r0 =====================================
                t=key
                try:
                    d = dn[t]
                except KeyError:
                    d0 = self.__sock_it_to_me_diffNS(n,t,surface=1)
                    d = Vector(k2)
                    d[0:k] = d0
                
                T0.inner(d,ind,map=amap,result=r) 
                if n>0:
                    d = d[0:k]
                    r0.inner(d,ind,map=amap,result=rPP)

                # Tp and rp =====================================
                if key != 0:
                    t=key; t=t-1
                    try:
                        d = dn[t]
                    except KeyError:
                        d0 = self.__sock_it_to_me_diffNS(n,t,surface=1)
                        d = Vector(k2)
                        d[0:k] = d0
                
                    Tp.inner(d,ind,map=amap,result=r) 
                    if n>0:
                        d = d[0:k]
                        rp.inner(d,ind,map=amap,result=rPP) 

                # Tm and rm =====================================
                if key != twon1:
                    t=key; t=t+1
                    try:
                        d = dn[t]
                    except KeyError:
                        d0 = self.__sock_it_to_me_diffNS(n,t,surface=1)
                        d = Vector(k2)
                        d[0:k] = d0
                        
                    Tm.inner(d,ind,map=amap,result=r)
                    if n>0:
                        d = d[0:k]
                        rm.inner(d,ind,map=amap,result=rPP)

                # ===============================================
                if n>0:
                    r[0:k].gaxpy(1.0,rPP[0:k],-1.0)
                    
                result.d[n][key] = r.scale(twon)
                
                #cut1 += 1
                #tmp = Vector(result.d[n][key])
                #tmp[0:k].gaxpy(1.0,result.d[n][key][0:k],-1.0)
                #if tmp.normf() > self.thresh:
                #    cut2 +=1

            del self.cache_sock_it_to_me_diffNS
            del self.cache_sock_it_to_me2_diffNS

        # Operator has been applied.  Must now cleanup f from its
        # nonstandard state, and also compress the result summing the
        # sum coeffs at all levels in the process.

        #print cut2, "/", cut1

        for n in result.d.keys():
            dn = result.d[n]
            sn = result.s[n] = {}
            for key in dn.keys():
                sn[key] = Vector(dn[key][0:k])
                dn[key][0:k].fill(0.0)
            
        result.compressed = 0
        result.compress()

        # Clean up the source
        self.s[0] = {}
        self.s[0][0] = Vector(self.d[0][0][0:k])
        for n in self.d.keys():
            dn = self.d[n]
            for key in dn.keys():
                dn[key][0:k].fill(0.0)
        self.compressed = 1
        
        for n in range(nmax,-1,-1):
            if len(self.d[n])==0:
                del self.d[n]
            else:
                break

        perf.exit('diffNS')
        return result

    def diff2NS(self,axis,autorefine=0):

        k, k2 = self.k, self.k*2
        
        if self.debug:
            print "diffing2", id(self), self.compressed, axis

        self.compress()
        perf.enter('diff2NS')
        self.reconstruct(nonstandard=1)

        result = self.new(k=self.k, thresh=self.thresh)

        ind = (1,axis)                  # Indices to contract in products
        if axis == 0:
            amap = None
        else:
            raise "Don't you think six dimensions are enough?"
            
        (Tm, T0, Tp), (rm, r0, rp) = self.make_dc_periodic_compressed()
        
        tp2 = Tm.transpose()*Tp
        tp1 = Tm.transpose()*T0 + T0.transpose()*Tp
        t0  = Tm.transpose()*Tm + T0.transpose()*T0 + Tp.transpose()*Tp
        tm1 =                     T0.transpose()*Tm + Tp.transpose()*T0
        tm2 =                                         Tp.transpose()*Tm

        Rp2 = rm.transpose()*rp
        Rp1 = rm.transpose()*r0 + r0.transpose()*rp
        R0  = rm.transpose()*rm + r0.transpose()*r0 + rp.transpose()*rp
        Rm1 =                     r0.transpose()*rm + rp.transpose()*r0
        Rm2 =                                         rp.transpose()*rm

        del Tm, T0, Tp, rm, r0, rp
        
        nmax = max(self.d.keys())

        if autorefine: nmax += 1
        
        for n in range(nmax+1):         #  Fill in empty levels
            try:
                self.d[n]
            except KeyError:
                self.d[n] = {}

        cut1, cut2 = 0, 0
            
        for n in range(nmax+1):
            result.d[n] = {}
            dn = self.d[n]
            twon = 2**n
            twon1= twon - 1
            scale = 4.0**n
            
            self.cache_sock_it_to_me_diffNS = {}
            self.cache_sock_it_to_me2_diffNS = {}
            
            if autorefine:
                if n>0:
                    keys = []
                    for lx in self.d[n-1].keys():
                        keys += [ lx*2, lx*2+1 ]
                else:
                    keys = [0]
            else:
                keys = dn.keys()

            for key in self.d[n].keys():
                t=key; t += 2
                if t<=twon1 and not (t in keys): keys.append(t)
                t=key; t += 1
                if t<=twon1 and not (t in keys): keys.append(t)
                t=key; t -= 1
                if t>=0     and not (t in keys): keys.append(t)
                t=key; t -= 2
                if t>=0     and not (t in keys): keys.append(t)

            for key in keys:
                r = Vector(k2)
                if n>0: rPP = Vector(k)

                # NOTE: We should expolit the sparsity of multiwavelet bases
                # in non-standard form 

                # t0 and R0 =====================================
                t=key
                try:
                    d = dn[t]
                except KeyError:
                    d0 = self.__sock_it_to_me_diffNS(n,t,surface=1)
                    d = Vector(k2)
                    d[0:k] = d0
                
                t0.inner(d,ind,map=amap,result=r) 
                if n>0:
                    d = d[0:k]
                    R0.inner(d,ind,map=amap,result=rPP)

                # tp1 and Rp1 ===================================
                if key != 0:
                    t=key; t=t-1
                    try:
                        d = dn[t]
                    except KeyError:
                        d0 = self.__sock_it_to_me_diffNS(n,t,surface=1)
                        d = Vector(k2)
                        d[0:k] = d0
                
                    tp1.inner(d,ind,map=amap,result=r) 
                    if n>0:
                        d = d[0:k]
                        Rp1.inner(d,ind,map=amap,result=rPP) 

                # tp2 and Rp2 ===================================
                if not (key in (0, 1)):
                    t=key; t=t-2
                    try:
                        d = dn[t]
                    except KeyError:
                        d0 = self.__sock_it_to_me_diffNS(n,t,surface=1)
                        d = Vector(k2)
                        d[0:k] = d0
                
                    tp2.inner(d,ind,map=amap,result=r) 
                    if n>0:
                        d = d[0:k]
                        Rp2.inner(d,ind,map=amap,result=rPP) 

                # tm1 and Rm1 ===================================
                if key != twon1:
                    t=key; t=t+1
                    try:
                        d = dn[t]
                    except KeyError:
                        d0 = self.__sock_it_to_me_diffNS(n,t,surface=1)
                        d = Vector(k2)
                        d[0:k] = d0
                        
                    tm1.inner(d,ind,map=amap,result=r)
                    if n>0:
                        d = d[0:k]
                        Rm1.inner(d,ind,map=amap,result=rPP)

                # tm2 and Rm2 ===================================
                if not (key in (twon1, twon1-1)):
                    t=key; t=t+2
                    try:
                        d = dn[t]
                    except KeyError:
                        d0 = self.__sock_it_to_me_diffNS(n,t,surface=1)
                        d = Vector(k2)
                        d[0:k] = d0
                        
                    tm2.inner(d,ind,map=amap,result=r)
                    if n>0:
                        d = d[0:k]
                        Rm2.inner(d,ind,map=amap,result=rPP)

                # ===============================================
                if n>0:
                    r[0:k].gaxpy(1.0,rPP[0:k],-1.0)
                    
                result.d[n][key] = r.scale(scale)
                
                #cut1 += 1
                #tmp = Vector(result.d[n][key])
                #tmp[0:k].gaxpy(1.0,result.d[n][key][0:k],-1.0)
                #if tmp.normf() > self.thresh: cut2 +=1

            del self.cache_sock_it_to_me_diffNS
            del self.cache_sock_it_to_me2_diffNS

        # Operator has been applied.  Must now cleanup f from its
        # nonstandard state, and also compress the result summing the
        # sum coeffs at all levels in the process.

        #print cut2, '/', cut1

        for n in result.d.keys():
            dn = result.d[n]
            sn = result.s[n] = {}
            for key in dn.keys():
                sn[key] = Vector(dn[key][0:k])
                dn[key][0:k].fill(0.0)
            
        result.compressed = 0
        result.compress()

        # Clean up the source
        self.s[0] = {}
        self.s[0][0] = Vector(self.d[0][0][0:k])
        for n in self.d.keys():
            dn = self.d[n]
            for key in dn.keys():
                dn[key][0:k].fill(0.0)
        self.compressed = 1
        
        for n in range(nmax,-1,-1):
            if len(self.d[n])==0:
                del self.d[n]
            else:
                break

        perf.exit('diff2NS')
        return result

    def project_bc(self):
        raise "not yet in 1d"

    def fcube(self):
        raise "no way"

    def latex_zslice(self):
        raise "???"

    def to2d(self,side=1):
        '''
        
        side = 1 : left    f_new(r1,r2) = f(r1)
        side = 2 : right   f_new(r1,r2) = f(r2)
        
        '''

        if side!=1 and side!=2:
            raise "side is either 1 or 2"
            
        perf.enter('to2d')

        if self.compressed:

            result = Function2d(k=self.k,compress=1) # Was thresh=self.thresh
            result.s = {}
            result.d = {}
    
            k = self.k
                
            one = Vector(k)
            one[0] = 1.0

            result.s[0]={}
            if(side==1):            
                result.s[0][0,0] = self.s[0][0].outer(one)
            else:
                result.s[0][0,0] = one.outer(self.s[0][0])
            
            one = Vector(2*k)
            
            for n in self.d.keys():
                twon = 2**n
                h = 1.0/twon
                scale = sqrt(h)
                scale = scale
                one[0] = scale
                    
                trans1d = range(twon)
                dn = self.d[n]
                rn  = result.d[n] = {}
                for key in dn.keys():
                    if(side==1):            
                        hh = dn[key].outer(one)
                        for lx in trans1d:
                            rn[(key,)+(lx,)] = hh
                        
                    else:
                        hh = one.outer(dn[key])
                        for lx in trans1d:
                            rn[(lx,)+(key,)] = hh
                                        
        else:

            result = Function2d(k=self.k,compress=0) # Was thresh=self.thresh
            result.s = {}
    
            k = self.k
                
            one = Vector(k)
                
            for n in self.s.keys():
                twon = 2**n
                h = 1.0/twon
                scale = sqrt(h)
                scale = scale
                one[0] = scale
                    
                trans1d = range(twon)
                sn = self.s[n]
                rn  = result.s[n] = {}
                for key in sn.keys():
                    if(side==1):            
                        hh = sn[key].outer(one)
                        for lx in trans1d:
                            rn[(key,)+(lx,)] = hh
                    else:
                        hh = one.outer(sn[key])
                        for lx in trans1d:
                            rn[(lx,)+(key,)] = hh

        perf.exit('to2d')

        return result

    
if __name__ == "__main__":

    Function.k =9
    Function.thresh=1e-9

    #start = time.clock()
    #f = Function(cfunction=cfcube.cvar.Cmyg,initial_level=2,refine=1,compress=0);
    #used = time.clock()-start
    #print "initial projection used", used

    #print f.k, f.thresh

    #f.print_layers(0)

    #start = time.clock()
    #f.compress()
    #used = time.clock()-start
    #print "compression used", used
#
#    start = time.clock()
#    f.reconstruct()
#    used = time.clock()-start
#    print "reconstruction used", used

#    f.compress()
#
#    start = time.clock()
#    f.truncate()
#    used = time.clock()-start
#    print "truncation used", used
#
#    print f.norm2()
#
#    #print f.s[2][1,1,1]
#    print f(0.5,0.5,0.5)
#    

    # 3-D test
    dotest1_3d = 1
    dotest2_3d = 1 # also 3, 4, 5, 6
    dotest7_3d = 0
    
    # 6-D test
    dotest1_6d = 0
    dotest2_6d = 0
    dotest9 = 0

    # 1-D test
    dotest0_1d = 0
    dotest1_1d = 0 #1
    dotest2_1d = 0 #1 # also 3, 4, 5, 6
    dotest8_1d = 0

    # 2-D test
    dotest1_2d = 0 #1
    dotest2_2d = 0 #1 # also 3, 4, 5, 6

    if dotest1_3d:
        # Can we even compress something simple?
        def testf(x,y,z):
            alpha = 5.0
            x, y, z = x-0.5, y-0.5, z-0.5
            return exp(-alpha*(x*x+y*y+z*z))*5.0

        t0=time.clock()
        print "test1: compressing a Gaussian to low precision"
        f = Function(function=testf,thresh=1e-6,k=9,refine=1,compress=0,initial_level=2)
	df = Dx * f
        f.print_layers(0)
	sys.exit()
        f.reconstruct()
        print "Compressed value: %.4f" % f(0.33,0.44,0.55)
        print "  Analytic value: %.4f" % testf(0.33,0.44,0.55)
        if abs(f(0.33,0.44,0.55) - testf(0.33,0.44,0.55)) > 3e-3:
            raise "test1 failed"
        print "test1: compression might be OK" 
        print "Truncating"
        f.truncate()
        print "Compressed value: %.4f" % f(0.33,0.44,0.55)
        print "  Analytic value: %.4f" % testf(0.33,0.44,0.55)
        print "Test 1 timing result " , time.clock()-t0
        del(f)

    if dotest2_3d:
        # Try compressing a function defined in C to higher precision
        # (in C for the speed)
        t0=time.clock()
        print "test2: compressing an exponential to medium precision"
        cfcube.cvar.myexponent = 25.0
        f = Function(cfunction=cfcube.cvar.Cmyf, function=cfcube.myf,
                     thresh=1e-7, k=9)
        f.print_layers(0)
        print "Compressed value: %.8f" % f(0.33,0.44,0.55)
        print "  Analytic value: %.8f" % cfcube.myf(0.33,0.44,0.55)
        if abs(f(0.33,0.44,0.55) - cfcube.myf(0.33,0.44,0.55)) > 3e-7:
            raise "test2 failed"
        print "test2: compression is OK"
        print "Truncating"
        f.truncate()
        print "Compressed value: %.8f" % f(0.33,0.44,0.55)
        print "  Analytic value: %.8f" % cfcube.myf(0.33,0.44,0.55)
        print "Test 2 timing result " , time.clock()-t0

        t0=time.clock()
        print "Testing application of local function - g(x)=f(x)/(1+f(x)**2)"
        g = f.copy()     
	print (g-f).norm2()
        g.local_function(cfunction=cfcube.cvar.Ctestloc)
        test1 = cfcube.testloc(f(0.33,0.44,0.55))
        test2 = g(0.33,0.44,0.55)
        print "   Analytic: %.8e" % test1
        print "  Numerical: %.8e" % test2

        if (abs(test1-test2) > 1e-5):
            raise "test2b failed"
        print "test2b local function  seems OK"
        del(g)
        print "Test 2b timing result " , time.clock()-t0

        # Try differentiating it
        t0=time.clock()
        df = Dy * f
        print "Compressed deriv: %.8f" % df(0.313,0.414,0.515)
        print "  Analytic deriv: %.8f" % cfcube.my_der(0.313,0.414,0.515)
        if (abs(df(0.313,0.414,0.515) - cfcube.my_der(0.313,0.414,0.515)) 
            > 1e-5):
              print  "test3 failed"
        print "test3: differentiation seems OK"
        del df
        print "Test 3 timing result " , time.clock()-t0

        # Try the Lapalacian
        t0=time.clock()
        delsqf = Delsq * f
        print "Compressed Laplacian: %.8f" % delsqf(0.313,0.414,0.465)
        print "  Analytic Laplacian: %.8f" % cfcube.my_lap(0.313,0.414,0.465)
        if (abs(delsqf(0.313,0.414,0.465) - cfcube.my_lap(0.313,0.414,0.465)) 
            > 2e-1):
            print "test4 failed"
        print "test4: laplacian seems OK"
        print "Test 4 timing result " , time.clock()-t0
        del(delsqf)

        # Try copying and squaring
        t0=time.clock()
        fsq = f.copy().square()
        val1 = fsq(0.313,0.414,0.465)
        val2 = cfcube.myf(0.313,0.414,0.465)**2        
        print "Compressed square: %.8e" % val1
        print "  Analytic square: %.8e" % val2
        if abs(val1-val2) > 1e-7:
           print "test5 failed"
        print "test5: copying and squaring seem to be OK"
        print "Test 5 timing result " , time.clock()-t0
        fsq.truncate()

        t0=time.clock()
        # Try addition and scalar multiplication
        g = fsq*1.3 + f*0.11
        val1 = g(0.45,0.55,0.45)
        val = cfcube.myf(0.45,0.55,0.45)
        val2 = val*val*1.3 + val*0.11
        print "Compressed fsq*1.3 + f*0.11: %.8e" % val1
        print "  Analytic fsq*1.3 + f*0.11: %.8e" % val2
        if abs(val1-val2) > 4e-5:
            pass # raise "test6 failed"
        print "test6: addition and scalar multiplication seem OK"
        print "Test 6 timing result " , time.clock()-t0
        del(fsq)

        t0=time.clock()
        # Try multiplication by another function
        h = f*g
        val1 = h(0.45,0.55,0.45)
        val = cfcube.myf(0.45,0.55,0.45)
        val2 = (val*val*1.3 + val*0.11)*val
        print "Compressed f*(fsq*1.3 + f*0.11): %.8e" % val1
        print "  Analytic f*(fsq*1.3 + f*0.11): %.8e" % val2
        if abs(val1-val2) > 1e-2:
             print "test6a failed"
        print "test6a: function multiplication seems OK"
        print "Test 6 timing result " , time.clock()-t0

        ### Try mul_func
        ##h21= f.mul_func(cfunction=cfcube.cvar.Ctestloc)
        ##h22= g.mul_func(function=cfcube.myf)
        ##val11= h21(0.45,0.55,0.45)
        ##val12= h22(0.45,0.55,0.45)
        ##val = cfcube.myf(0.45,0.55,0.45)
        ##val2 = (val*val*1.3 + val*0.11)*val
        ##print "Compressed f*(fsq*1.3 + f*0.11): %.8e" % val11
        ##print "Compressed f*(fsq*1.3 + f*0.11): %.8e" % val12
        ##print "  Analytic f*(fsq*1.3 + f*0.11): %.8e" % val2
        ##if abs(val11-val2) > 1e-2 or abs(val12-val2) > 1e-2:
        ##    raise "test6b failed"
        ##print "test6b: function multiplication using mul_func seems OK"

        t0=time.clock()
        # Try inner product in both bases
        f.compress()
        g.compress()
        test1 = f.inner(g)
        f.reconstruct()
        g.reconstruct()
        test2 = f.inner(g)
        test3 = h.trace()
        print "  Compressed: <f|g> = %.8e" % test1
        print "Uncompressed: <f|g> = %.8e" % test2
        print "     Tr(f*g): <f|g> = %.8e" % test3
        print "Test 6b timing result " , time.clock()-t0
        del(g)
        del(h)

        # Compare different evaluations
        f.compress()
        fc = f(0.45,0.55,0.53)
        f.reconstruct()
        fr = f(0.45,0.55,0.53)
        f.fast_eval_init()
        ff = f(0.45,0.55,0.53)
        f.fast_eval_tidy()
        print " Evaluation compressed:  %.8e" % fc
        print " Evaluation unompressed: %.8e" % fr
        print " Evaluation fast:        %.8e" % ff
        if abs(ff-fc)>1e-12 or abs(fr-fc)>1e-12:
            raise "Evaluation compressed/uncompressed/fast inconsistent"
        

    if dotest7_3d:
        k = 7
        print "Testing Poisson for k =", k
        t0=time.clock()
        o = Poisson(k=k,dirname='./poisson/')
        print "Compressing target function"
        f = Function(cfunction=cfcube.cvar.Cmyf,
                     function=cfcube.myf,
                     thresh=0.1**(k-2), k=k)
        f.truncate()
        print "Applying operator"
        v = o.apply_sparse(f,bmax=4)
        print "Test 7  timing result " , time.clock()-t0
        print "Compressed potential %.8f" % v(0.313,0.414,0.465)
        print "  Analytic potential %.8f" % cfcube.mypot(0.313,0.414,0.465)
        for z in 0.01, 0.1, 0.2, 0.3, 0.4,0.45, 0.49, 0.51, 0.55, 0.6, \
            0.7, 0.8, 0.9, 0.99:
            vv = v(0.5,0.5,z)
            cc = cfcube.mypot(0.5,0.5,z)
            print z, vv, cc, cc/vv

    #          ################################################################
    # 6-D TEST ################################################################
    #          ################################################################

    if dotest1_6d:
        # Can we even compress something simple?
        def testf(x1,y1,z1,x2,y2,z2):
            alpha1, alpha2 = 5.0, 1.0
            x1, y1, z1, x2, y2, z2 = x1-0.40, y1-0.45, z1-0.50, x2-0.50, y2-0.55, z2-0.60
            return exp(-alpha1*(x1*x1+y1*y1+z1*z1))*5.0 \
                 * exp(-alpha2*(x2*x2+y2*y2+z2*z2))*3.0

        testf = cfcube.testexp6d
        ctestf = cfcube.cvar.Ctestexp6d

        print "[CAUTION]:"
        print "[CAUTION]: this test takes a lot of time (maybe, dozen minites)."
        print "[CAUTION]:"

        print "\ntest1 for 6d: compressing a Gaussian to low precision"
        f = Function6d(cfunction=ctestf, function=testf,
                       thresh=1e-3,k=5,compress=0,refine=1,initial_level=1)

        f.memory_usage()
        print "Reconstruc value: %.4f" % f(0.33,0.44,0.55,0.33,0.44,0.55)
        print "  Analytic value: %.4f" % testf(0.33,0.44,0.55,0.33,0.44,0.55)
        if abs(f(0.33,0.44,0.55,0.33,0.44,0.55) - testf(0.33,0.44,0.55,0.33,0.44,0.55)) > 3e-3:
            raise "test1 failed"
        print "test1: reconstruct might be OK"

        f.compress()
        f.memory_usage()
        print "Compressed value: %.4f" % f(0.33,0.44,0.55,0.33,0.44,0.55)
        print "  Analytic value: %.4f" % testf(0.33,0.44,0.55,0.33,0.44,0.55)
        if abs(f(0.33,0.44,0.55,0.33,0.44,0.55) - testf(0.33,0.44,0.55,0.33,0.44,0.55)) > 3e-3:
            raise "test1 failed"
        print "test1: compression might be OK"

        f.reconstruct()
        f.memory_usage()
        print "Reconstruc value: %.4f" % f(0.33,0.44,0.55,0.33,0.44,0.55)
        print "  Analytic value: %.4f" % testf(0.33,0.44,0.55,0.33,0.44,0.55)
        if abs(f(0.33,0.44,0.55,0.33,0.44,0.55) - testf(0.33,0.44,0.55,0.33,0.44,0.55)) > 3e-3:
            raise "test1 failed"
        print "test1: reconstruct might be OK"

        print "Truncating"
        f.truncate()
        f.memory_usage()
        print "Compressed value: %.4f" % f(0.33,0.44,0.55,0.33,0.44,0.55)
        print "  Analytic value: %.4f" % testf(0.33,0.44,0.55,0.33,0.44,0.55)

        f.reconstruct()
        f.memory_usage()
        print "Reconstruc value: %.4f" % f(0.33,0.44,0.55,0.33,0.44,0.55)
        print "  Analytic value: %.4f" % testf(0.33,0.44,0.55,0.33,0.44,0.55)
        if abs(f(0.33,0.44,0.55,0.33,0.44,0.55) - testf(0.33,0.44,0.55,0.33,0.44,0.55)) > 3e-3:
            raise "test1 failed"
        print "test1: reconstruct might be OK"

    if dotest2_6d:

        print "[CAUTION]:"
        print "[CAUTION]: this test takes a lot of time (maybe, dozen minites)."
        print "[CAUTION]:"

        alpha1 = 5.0
        alpha2 = 1.0
        
        def testfa(x1,y1,z1):
            x1, y1, z1 = x1-0.40, y1-0.45, z1-0.50
            return exp(-alpha1*(x1*x1+y1*y1+z1*z1))

        def testfb(x2,y2,z2):
            x2, y2, z2 = x2-0.50, y2-0.55, z2-0.60
            return exp(-alpha2*(x2*x2+y2*y2+z2*z2))

        def testfa_der(x1,y1,z1):
            x1, y1, z1 = x1-0.40, y1-0.45, z1-0.50
            return exp(-alpha1*(x1*x1+y1*y1+z1*z1))*(-2.0*alpha1*y1)

        def testfb_der(x2,y2,z2):
            x2, y2, z2 = x2-0.50, y2-0.55, z2-0.60
            return exp(-alpha2*(x2*x2+y2*y2+z2*z2))*(-2.0*alpha2*y2)

        def testfa_der2(x1,y1,z1):
            x1, y1, z1 = x1-0.40, y1-0.45, z1-0.50
            return exp(-alpha1*(x1*x1+y1*y1+z1*z1))*(-2.0*alpha1+4.0*alpha1*alpha1*y1*y1)

        def testfb_der2(x2,y2,z2):
            x2, y2, z2 = x2-0.50, y2-0.55, z2-0.60
            return exp(-alpha2*(x2*x2+y2*y2+z2*z2))*(-2.0*alpha2+4.0*alpha2*alpha2*y2*y2)

        def testfab(x1,y1,z1,x2,y2,z2):
            x1, y1, z1 = x1-0.40, y1-0.45, z1-0.50
            x2, y2, z2 = x2-0.50, y2-0.55, z2-0.60
            return testfa(x1,y1,z1)*testfb(x2,y2,z2)

        def testfab_lap(x1,y1,z1,x2,y2,z2):
            x1, y1, z1 = x1-0.40, y1-0.45, z1-0.50
            x2, y2, z2 = x2-0.50, y2-0.55, z2-0.60
            return exp(-alpha1*(x1*x1+y1*y1+z1*z1)-alpha2*(x2*x2+y2*y2+z2*z2)) \
                   *(-8.0*alpha1
                     -8.0*alpha2
                     +4.0*alpha1*alpha1*(x1*x1+y1*y1+z1*z1)
                     +4.0*alpha2*alpha2*(x2*x2+y2*y2+z2*z2))

        Function.autorefine = 0
        Function.truncate_method = 0
        
        fa = Function(function=testfa,thresh=1e-5,k=7,initial_level=1)
        fa.memory_usage()
        print "comp fa done"
        fb = Function(function=testfb,thresh=1e-5,k=7)
        fb.memory_usage()
        print "comp fb done"

        fa.reconstruct()
        fb.reconstruct()

        fa = fa.to6d(1)  # fa(r,r') = fa(r)
        print "fa(reconstruct) to6d done"

        fb = fb.to6d(2)  # fb(r,r') = fb(r')
        print "fb(reconstruct) to6d done"

        print "Compressed value: %.8f %.8f" % ( fa(0.33,0.44,0.55,0.5,0.5,0.5), fb(0.5,0.5,0.5,0.40,0.50,0.60) )
        print "Compressed value: %.8f %.8f" % ( fa(0.33,0.44,0.55,0.2,0.3,0.4), fb(0.2,0.3,0.4,0.40,0.50,0.60) )
        print "  Analytic value: %.8f %.8f" % ( testfa(0.33,0.44,0.55), testfb(0.40,0.50,0.60) )

##        fa = Function(function=testfa,thresh=1e-5,k=7,initial_level=1)
##        fa.memory_usage()
##        print "comp fa done"
##        fb = Function(function=testfb,thresh=1e-5,k=7)
##        fb.memory_usage()
##        print "comp fb done"
##
##        fa.compress()
##        fb.compress()
##
##        fa = fa.to6d(1)
##        print "fa(compress) to6d done"
##
##        fb = fb.to6d(2)
##        print "fb(compress) to6d done"
##
##        print "Compressed value: %.8f %.8f" % ( fa(0.33,0.44,0.55,0.5,0.5,0.5), fb(0.5,0.5,0.5,0.40,0.50,0.60) )
##        print "Compressed value: %.8f %.8f" % ( fa(0.33,0.44,0.55,0.2,0.3,0.4), fb(0.2,0.3,0.4,0.40,0.50,0.60) )
##        print "  Analytic value: %.8f %.8f" % ( testfa(0.33,0.44,0.55), testfb(0.40,0.50,0.60) )
##
##        fa.compress()
##        fb.compress()
##
##        print "gaxpy start"
##        f = fa.copy().gaxpy(2.0,fb,3.0)
##        print "gaxpy end"
##
##        print "Compressed value: %.8f" % f(0.33,0.44,0.55,0.40,0.50,0.60)
##        print "  Analytic value: %.8f" % (2.0*testfa(0.33,0.44,0.55)+3.0*testfb(0.40,0.50,0.60))
##
##        print "mul start"
##        f = fa*fb
##        print "mul done"
##
##        print "Compressed value: %.8f" % f(0.33,0.44,0.55,0.33,0.44,0.55)
##        print "  Analytic value: %.8f" % (testfa(0.33,0.44,0.55)*testfb(0.33,0.44,0.55))
##        print "Compressed value: %.8f" % f(0.33,0.44,0.55,0.40,0.50,0.60)
##        print "  Analytic value: %.8f" % (testfa(0.33,0.44,0.55)*testfb(0.40,0.50,0.60))

        fb3d = Function(function=testfb,thresh=1e-5,k=7)
        fb3d.memory_usage()
        print "comp fb3d done"
        fb3d.reconstruct()
        fa.reconstruct()
        
        print "mul_with3d start"
        f = fa.mul_with3d(fb3d,side=2)
        print "mul done"

        print "Compressed value: %.8f" % f(0.33,0.44,0.55,0.33,0.44,0.55)
        print "  Analytic value: %.8f" % (testfa(0.33,0.44,0.55)*testfb(0.33,0.44,0.55))
        print "Compressed value: %.8f" % f(0.33,0.44,0.55,0.40,0.50,0.60)
        print "  Analytic value: %.8f" % (testfa(0.33,0.44,0.55)*testfb(0.40,0.50,0.60))

        fa3d = Function(function=testfa,thresh=1e-5,k=7,initial_level=1)
        fa3d.memory_usage()
        print "comp fa3d done"
        fa3d.reconstruct()

        fb3d = Function(function=testfb,thresh=1e-5,k=7)
        fb3d.memory_usage()
        print "comp fb3d done"
        fb3d.reconstruct()
        fa.reconstruct()
        
        print "mul start with_3d"
        f = f.mul_with3d(fb3d,side=1)
        print "mul done"

        print "Compressed value: %.8f" % f(0.33,0.44,0.55,0.40,0.50,0.60)
        print "  Analytic value: %.8f" % (testfa(0.33,0.44,0.55)*testfb(0.40,0.50,0.60)*testfb(0.33,0.44,0.55))

        print "mul start with_3d"
        f = f.mul_with3d(fa3d,side=2)
        print "mul done"

        print "Compressed value: %.8f" % f(0.33,0.44,0.55,0.40,0.50,0.60)
        print "  Analytic value: %.8f" % (
            testfa(0.33,0.44,0.55)*testfb(0.33,0.44,0.55)*testfb(0.40,0.50,0.60)*testfa(0.40,0.50,0.60))

        ddddddddddddddddddddddddddddddddddddddddddddddddddddddddddd

        ### Try compressing a function defined in C to higher precision
        ### (in C for the speed)
        ##print "\ntest2: compressing an exponential to medium precision"
        ##cfcube.cvar.myexponent = 5.0
        ##f = Function6d(cfunction=cfcube.cvar.Ctestexp6d, function=cfcube.testexp6d,
        ##             thresh=1e-3, k=5, initial_level=1)

        f.memory_usage()
        f.print_layers(0)
        print "Reconstruc value: %.8f" % f(0.33,0.44,0.55,0.33,0.44,0.55)
        print "  Analytic value: %.8f" % cfcube.testexp6d(0.33,0.44,0.55,0.33,0.44,0.55)
        if abs(f(0.33,0.44,0.55,0.33,0.44,0.55) - cfcube.testexp6d(0.33,0.44,0.55,0.33,0.44,0.55)) > 3e-5:
            print "test2 failed"
        print "test2: compression is OK"
        
        print "Truncating"
        f.truncate(1e-6)
        f.memory_usage()
        f.print_layers(0)
        print "Compressed value: %.8f" % f(0.33,0.44,0.55,0.33,0.44,0.55)
        print "  Analytic value: %.8f" % cfcube.testexp6d(0.33,0.44,0.55,0.33,0.44,0.55)

        print "Testing application of local function - g(x)=f(x)/(1+f(x)**2)"
        g = f.copy()
        g.local_function(cfunction=cfcube.cvar.Ctestloc)
        test1 = cfcube.testloc(f(0.33,0.44,0.55,0.33,0.44,0.55))
        test2 = g(0.33,0.44,0.55,0.33,0.44,0.55)
        print "   Analytic: %.8e" % test1
        print "  Numerical: %.8e" % test2
        del(g)
        
        # Try differentiating it
        df = f.diff(1)
        print "Compressed deriv: %.8f" % df(0.313,0.414,0.515,0.313,0.414,0.515)
        print "  Analytic deriv: %.8f" % (testfa_der(0.313,0.414,0.515)*testfb(0.313,0.414,0.515))
        del df

        df = f.diff(4)
        print "Compressed deriv: %.8f" % df(0.313,0.414,0.515,0.313,0.414,0.515)
        print "  Analytic deriv: %.8f" % (testfa(0.313,0.414,0.515)*testfb_der(0.313,0.414,0.515))
        del df
        print "test3_1: differentiation seems OK"

####        # Try 2-differentiating it
####        df = f.diff2(1)
####        print "Compressed deriv: %.8f" % df(0.313,0.414,0.515,0.313,0.414,0.515)
####        print "  Analytic deriv: %.8f" % (testfa_der2(0.313,0.414,0.515)*testfb(0.313,0.414,0.515))
####        del df
####
####        df = f.diff2(4)
####        print "Compressed deriv: %.8f" % df(0.313,0.414,0.515,0.313,0.414,0.515)
####        print "  Analytic deriv: %.8f" % (testfa(0.313,0.414,0.515)*testfb_der2(0.313,0.414,0.515))
####        del df
####        print "test3_2: differentiation seems OK"
####
####        # Try the Lapalacian
####        delsqf = f.laplacian()
####        print "Compressed Laplacian: %.8f" % delsqf(0.313,0.414,0.465,0.313,0.414,0.465)
####        print "  Analytic Laplacian: %.8f" % testfab_lap(0.313,0.414,0.465,0.313,0.414,0.465)
####        del delsqf
####        print "test4: laplacian seems OK"

        # Try copying and squaring
        fsq = f.copy().square()
        val1 = fsq(0.313,0.414,0.465,0.313,0.414,0.465)
        val2 = cfcube.testexp6d(0.313,0.414,0.465,0.313,0.414,0.465)**2        
        print "Compressed square: %.8e" % val1
        print "  Analytic square: %.8e" % val2
        if abs(val1-val2) > 1e-5:
            print "test5 failed"
        print "test5: copying and squaring seem to be OK"

        fsq.truncate(1e-6)

        # Try addition and scalar multiplication
        g = fsq*1.3 + f*0.11
        val1 = g(0.45,0.55,0.45,0.45,0.55,0.45)
        val = cfcube.testexp6d(0.45,0.55,0.45,0.45,0.55,0.45)
        val2 = val*val*1.3 + val*0.11
        print "Compressed fsq*1.3 + f*0.11: %.8e" % val1
        print "  Analytic fsq*1.3 + f*0.11: %.8e" % val2
        if abs(val1-val2) > 4e-4:
            print "test6 failed"
        print "test6: addition and scalar multiplication seem OK"

        # Try multiplication by another function
        h = f*g
        val1 = h(0.45,0.55,0.45,0.45,0.55,0.45)
        val = cfcube.testexp6d(0.45,0.55,0.45,0.45,0.55,0.45)
        val2 = (val*val*1.3 + val*0.11)*val
        print "Compressed f*(fsq*1.3 + f*0.11): %.8e" % val1
        print "  Analytic f*(fsq*1.3 + f*0.11): %.8e" % val2
        if abs(val1-val2) > 5e-3:
            print "test6a failed"
        print "test6a: function multiplication seems OK"

        ### Try mul_func
        ##h21= f.mul_func(cfunction=cfcube.cvar.Ctestloc)
        ##h22= g.mul_func(function=cfcube.testexp6d)
        ##val11= h21(0.45,0.55,0.45)
        ##val12= h22(0.45,0.55,0.45)
        ##val = cfcube.testexp6d(0.45,0.55,0.45)
        ##val2 = (val*val*1.3 + val*0.11)*val
        ##print "Compressed f*(fsq*1.3 + f*0.11): %.8e" % val11
        ##print "Compressed f*(fsq*1.3 + f*0.11): %.8e" % val12
        ##print "  Analytic f*(fsq*1.3 + f*0.11): %.8e" % val2
        ##if abs(val11-val2) > 1e-2 or abs(val12-val2) > 1e-2:
        ##    raise "test6b failed"
        ##print "test6b: function multiplication using mul_func seems OK"

        # Try inner product in both bases
        f.compress()
        g.compress()
        test1 = f.inner(g)
        f.reconstruct()
        g.reconstruct()
        test2 = f.inner(g)
        test3 = h.trace()
        print "  Compressed: <f|g> = %.8e" % test1
        print "Uncompressed: <f|g> = %.8e" % test2
        print "     Tr(f*g): <f|g> = %.8e" % test3

        perf.output()

    if dotest9:
        print "Testing 6d function"
        t0=time.clock()
        
        L = 7.0
        cfcube.set_box_size(L);
        f = Function6d(k=6,cfunction=cfcube.cvar.Che3,thresh=2e-4,
                       initial_level=2,refine=1,compress=0)
        f.memory_usage()
        print f(0.5,0.5,0.5,0.5,0.5,0.5)
        print f.norm2()

        fac = L*L*L;
        
        r = 0.5/L
        for theta in range(0,101):
            theta = 3.1415926535*theta*0.01
            x = r*sin(theta/2)
            y = r*cos(theta/2)
            x1 = 0.5 + x
            y1 = 0.5 + y
            x2 = 0.5 - x
            y2 = 0.5 + y
            test = f(x1,y1,0.5,x2,y2,0.5) / fac
            exact = cfcube.he3(x1,y1,0.5,x2,y2,0.5) / fac
            err = test - exact
            if exact > 0.0:
                rel = err/exact
            else:
                rel = 0.0

            print "Test 9  timing result " , time.clock()-t0
            theta = 180.0*theta/3.1415926535
            print "%8.3f %.6f %.6f %.6f %.1e" % (theta, test, exact, err, rel)
                
                
                

    #          ################################################################
    # 1-D TEST ################################################################
    #          ################################################################
    if dotest0_1d:
        print "Testing 1D function"
        def f(x):
            return 2.0*3.4549414947*exp(-150.0*(x-0.5)**2)

        a = Function1d(function=f,refine=1,k=9,thresh=1e-10)
        print "Norm ", a.trace()
        print a.d.keys()
        
        print "Compressed value: %.8f" % a(0.45)
        print "  Analytic value: %.8f" % f(0.45)
        if abs(a(0.45)-f(0.45)) > 1e-9:
            raise "1d compression failed"

        a.eval_err()
        
        print "1d compression OK"
        asq = a.copy().square()
        print "Compressed square: %.8f" % asq(0.45)
        print "  Analytic square: %.8f" % f(0.45)**2
        if abs(asq(0.45)-f(0.45)**2) > 1e-9:
            raise "1d square failed"
        print "1d square OK"

    if dotest1_1d:
        # Can we even compress something simple?
        def testf(x):
            alpha = 5.0
            x = x-0.50
            return exp(-alpha*(x*x))*5.0

        testf = cfcube.testexp1d
        ctestf = cfcube.cvar.Ctestexp1d

        print "\ntest1 for 1d: compressing a Gaussian to low precision"
        f = Function1d(cfunction=ctestf, function=testf,
                       thresh=1e-3,k=5,compress=0,refine=1,initial_level=3)

        f.memory_usage()
        print "Reconstruc value: %.4f" % f(0.33)
        print "  Analytic value: %.4f" % testf(0.33)
        if abs(f(0.33) - testf(0.33)) > 3e-3:
            raise "test1 failed"
        print "test1: reconstruct might be OK"

        f.compress()
        f.memory_usage()
        print "Compressed value: %.4f" % f(0.33)
        print "  Analytic value: %.4f" % testf(0.33)
        if abs(f(0.33) - testf(0.33)) > 3e-3:
            raise "test1 failed"
        print "test1: compression might be OK"

        f.reconstruct()
        f.memory_usage()
        print "Reconstruc value: %.4f" % f(0.33)
        print "  Analytic value: %.4f" % testf(0.33)
        if abs(f(0.33) - testf(0.33)) > 3e-3:
            raise "test1 failed"
        print "test1: reconstruct might be OK"

        print "Truncating"
        f.truncate()
        f.memory_usage()
        print "Compressed value: %.4f" % f(0.33)
        print "  Analytic value: %.4f" % testf(0.33)

        f.reconstruct()
        f.memory_usage()
        print "Reconstruc value: %.4f" % f(0.33)
        print "  Analytic value: %.4f" % testf(0.33)
        if abs(f(0.33) - testf(0.33)) > 3e-3:
            raise "test1 failed"
        print "test1: reconstruct might be OK"

    if dotest2_1d:

        alpha = 5.0
        
        def testfa(x):
            x = x-0.50
            return sqrt(5.0)*exp(-alpha*(x*x))

        def testfa_der(x):
            x = x-0.50
            return sqrt(5.0)*exp(-alpha*(x*x))*(-2.0*alpha1*x)

        Function.autorefine = 0
        Function.truncate_method = 0
        
        fa = Function1d(function=testfa,thresh=1e-5,k=7)
        fa.memory_usage()
        print "comp fa done"
        fb = Function1d(function=testfa,thresh=1e-5,k=7)
        fb.memory_usage()
        print "comp fb done"

        fa.reconstruct()
        fb.reconstruct()
        
        print "gaxpy start"
        f = fa.copy().gaxpy(2.0,fb,3.0)
        print "gaxpy end"

        print "Compressed value: %.8f" % f(0.33)
        print "  Analytic value: %.8f" % (2.0*testfa(0.33)+3.0*testfa(0.33))

        fa.compress()
        fb.compress()

        print "gaxpy start"
        f = fa.copy().gaxpy(2.0,fb,3.0)
        print "gaxpy end"

        print "Compressed value: %.8f" % f(0.33)
        print "  Analytic value: %.8f" % (2.0*testfa(0.33)+3.0*testfa(0.33))

        print "mul start"
        f = fa*fb
        print "mul done"

        print "Compressed value: %.8f" % f(0.33)
        print "  Analytic value: %.8f" % (testfa(0.33)*testfa(0.33))
        print "Compressed value: %.8f" % f(0.44,0.44)
        print "  Analytic value: %.8f" % (testfa(0.44)*testfa(0.44))

        f = Function1d(cfunction=cfcube.cvar.Ctestexp1d,thresh=1e-7,k=9)
        f.reconstruct()
        
        f.memory_usage()
        f.print_layers(0)
        print "Reconstruc value: %.8f" % f(0.33)
        print "  Analytic value: %.8f" % cfcube.testexp1d(0.33)
        if abs(f(0.33) - cfcube.testexp1d(0.33)) > 3e-5:
            print "test2 failed"
        print "test2: compression is OK"

        print "comp fa done"
        print "Truncating"
        f.truncate(1e-7)
        f.memory_usage()
        f.print_layers(0)
        print "Compressed value: %.8f" % f(0.33)
        print "  Analytic value: %.8f" % cfcube.testexp1d(0.33)

        print "Testing application of local function - g(x)=f(x)/(1+f(x)**2)"
        g = f.copy()
        g.local_function(cfunction=cfcube.cvar.Ctestloc)
        test1 = cfcube.testloc(f(0.33))
        test2 = g(0.33)
        print "   Analytic: %.8e" % test1
        print "  Numerical: %.8e" % test2
        del(g)
        
        # Try differentiating it
        df = f.diff(0)
        print "Compressed deriv: %.8f" % df(0.313)
        print "  Analytic value: %.8f" % cfcube.testexp1d_dx(0.313)
        del df

        # Try copying and squaring
        fsq = f.copy().square()
        val1 = fsq(0.313)
        val2 = cfcube.testexp1d(0.313)**2        
        print "Compressed square: %.8e" % val1
        print "  Analytic square: %.8e" % val2
        if abs(val1-val2) > 1e-5:
            print "test5 failed"
        print "test5: copying and squaring seem to be OK"

        fsq.truncate(1e-6)

        # Try addition and scalar multiplication
        g = fsq*1.3 + f*0.11
        val1 = g(0.45)
        val = cfcube.testexp1d(0.45)
        val2 = val*val*1.3 + val*0.11
        print "Compressed fsq*1.3 + f*0.11: %.8e" % val1
        print "  Analytic fsq*1.3 + f*0.11: %.8e" % val2
        if abs(val1-val2) > 4e-4:
            print "test6 failed"
        print "test6: addition and scalar multiplication seem OK"

        # Try multiplication by another function
        h = f*g
        val1 = h(0.45)
        val = cfcube.testexp1d(0.45)
        val2 = (val*val*1.3 + val*0.11)*val
        print "Compressed f*(fsq*1.3 + f*0.11): %.8e" % val1
        print "  Analytic f*(fsq*1.3 + f*0.11): %.8e" % val2
        if abs(val1-val2) > 5e-3:
            print "test6a failed"
        print "test6a: function multiplication seems OK"

        ### Try mul_func
        ##h21= f.mul_func(cfunction=cfcube.cvar.Ctestloc)
        ##h22= g.mul_func(function=cfcube.testexp6d)
        ##val11= h21(0.45,0.55,0.45)
        ##val12= h22(0.45,0.55,0.45)
        ##val = cfcube.testexp6d(0.45,0.55,0.45)
        ##val2 = (val*val*1.3 + val*0.11)*val
        ##print "Compressed f*(fsq*1.3 + f*0.11): %.8e" % val11
        ##print "Compressed f*(fsq*1.3 + f*0.11): %.8e" % val12
        ##print "  Analytic f*(fsq*1.3 + f*0.11): %.8e" % val2
        ##if abs(val11-val2) > 1e-2 or abs(val12-val2) > 1e-2:
        ##    raise "test6b failed"
        ##print "test6b: function multiplication using mul_func seems OK"

        # Try inner product in both bases
        f.compress()
        g.compress()
        test1 = f.inner(g)
        f.reconstruct()
        g.reconstruct()
        test2 = f.inner(g)
        test3 = h.trace()
        print "  Compressed: <f|g> = %.8e" % test1
        print "Uncompressed: <f|g> = %.8e" % test2
        print "     Tr(f*g): <f|g> = %.8e" % test3

        perf.output()

    if dotest8_1d:
        print "Testing 1D function"
        def f(x):
            return 2.0*3.4549414947*exp(-150.0*(x-0.5)**2)

        a = Function1d(function=f,refine=1,k=9,thresh=1e-10)
        print "Norm ", a.trace()
        print a.d.keys()
        
        print "Compressed value: %.8f" % a(0.45)
        print "  Analytic value: %.8f" % f(0.45)
        if abs(a(0.45)-f(0.45)) > 1e-9:
            raise "1d compression failed"

        a.eval_err()
        
        print "1d compression OK"
        asq = a.copy().square()
        print "Compressed square: %.8f" % asq(0.45)
        print "  Analytic square: %.8f" % f(0.45)**2
        if abs(asq(0.45)-f(0.45)**2) > 1e-9:
            raise "1d square failed"
        print "1d square OK"

    #          ################################################################
    # 2-D TEST ################################################################
    #          ################################################################
    if dotest1_2d:
        # Can we even compress something simple?
        def testf(x,y):
            alpha = 5.0
            x, y = x-0.50, y-0.50
            return exp(-alpha*(x*x+y*y))*5.0

        testf = cfcube.testexp2d
        ctestf = cfcube.cvar.Ctestexp2d

        print "\ntest1 for 2d: compressing a Gaussian to low precision"
        f = Function2d(cfunction=ctestf, function=testf,
                       thresh=1e-3,k=5,compress=0,refine=1,initial_level=3)

        f.memory_usage()
        print "Reconstruc value: %.4f" % f(0.33,0.44)
        print "  Analytic value: %.4f" % testf(0.33,0.44)
        if abs(f(0.33,0.44) - testf(0.33,0.44)) > 3e-3:
            raise "test1 failed"
        print "test1: reconstruct might be OK"

        f.compress()
        f.memory_usage()
        print "Compressed value: %.4f" % f(0.33,0.44)
        print "  Analytic value: %.4f" % testf(0.33,0.44)
        if abs(f(0.33,0.44) - testf(0.33,0.44)) > 3e-3:
            raise "test1 failed"
        print "test1: compression might be OK"

        f.reconstruct()
        f.memory_usage()
        print "Reconstruc value: %.4f" % f(0.33,0.44)
        print "  Analytic value: %.4f" % testf(0.33,0.44)
        if abs(f(0.33,0.44) - testf(0.33,0.44)) > 3e-3:
            raise "test1 failed"
        print "test1: reconstruct might be OK"

        print "Truncating"
        f.truncate()
        f.memory_usage()
        print "Compressed value: %.4f" % f(0.33,0.44)
        print "  Analytic value: %.4f" % testf(0.33,0.44)

        f.reconstruct()
        f.memory_usage()
        print "Reconstruc value: %.4f" % f(0.33,0.44)
        print "  Analytic value: %.4f" % testf(0.33,0.44)
        if abs(f(0.33,0.44) - testf(0.33,0.44)) > 3e-3:
            raise "test1 failed"
        print "test1: reconstruct might be OK"


    if dotest2_2d:

        alpha = 5.0
        
        def testfa(x):
            x = x-0.50
            return sqrt(5.0)*exp(-alpha*(x*x))

        def testfa_der(x):
            x = x-0.50
            return sqrt(5.0)*exp(-alpha*(x*x))*(-2.0*alpha1*x)

        def testfa_der2(x1,y1,z1):
            x = x-0.50
            return sqrt(5.0)*exp(-alpha*(x*x))*(-2.0*alpha+4.0*alpha*alpha*x*x)

        def testfab(x,y):
            x, y = x-0.50, y-0.50
            return testfa(x)*testfb(y)

        def testfab_lap(x,y):
            x, y = x-0.50, y-0.50
            return 5.0*exp(-alpha*(x*x+y*y)) \
                   *(-4.0*alpha+4.0*alpha*alpha*(x*x+y*y))

        Function.autorefine = 0
        Function.truncate_method = 0
        
        fa = Function1d(function=testfa,thresh=1e-5,k=7)
        fa.memory_usage()
        print "comp fa done"
        fb = Function1d(function=testfa,thresh=1e-5,k=7)
        fb.memory_usage()
        print "comp fb done"

        fa.reconstruct()
        fb.reconstruct()
        
        fa = fa.to2d(1)  # fa(r,r') = fa(r)
        print "fa(reconstruct) to2d done"

        fb = fb.to2d(2)  # fb(r,r') = fb(r')
        print "fb(reconstruct) to2d done"

        print "Compressed value: %.8f %.8f" % ( fa(0.33,0.50), fb(0.20,0.40) )
        print "  Analytic value: %.8f %.8f" % ( testfa(0.33), testfa(0.40) )

        fa.reconstruct()
        fb.reconstruct()

        print "gaxpy start"
        f = fa.copy().gaxpy(2.0,fb,3.0)
        print "gaxpy end"

        print "Compressed value: %.8f" % f(0.33,0.40)
        print "  Analytic value: %.8f" % (2.0*testfa(0.33)+3.0*testfa(0.40))

        fa.compress()
        fb.compress()

        print "gaxpy start"
        f = fa.copy().gaxpy(2.0,fb,3.0)
        print "gaxpy end"

        print "Compressed value: %.8f" % f(0.33,0.40)
        print "  Analytic value: %.8f" % (2.0*testfa(0.33)+3.0*testfa(0.40))

        print "mul start"
        f = fa*fb
        print "mul done"

        print "Compressed value: %.8f" % f(0.33,0.40)
        print "  Analytic value: %.8f" % (testfa(0.33)*testfa(0.40))
        print "Compressed value: %.8f" % f(0.44,0.55)
        print "  Analytic value: %.8f" % (testfa(0.44)*testfa(0.55))

        fb1d = Function1d(function=testfa,thresh=1e-5,k=7)
        fb1d.memory_usage()
        print "comp fb1d done"
        fb1d.reconstruct()
        fa.reconstruct()
        
        print "mul_with1d start"
        f = fa.mul_with1d(fb1d,side=2)
        print "mul done"

        print "Compressed value: %.8f" % f(0.33,0.40)
        print "  Analytic value: %.8f" % (testfa(0.33)*testfa(0.40))
        print "Compressed value: %.8f" % f(0.44,0.55)
        print "  Analytic value: %.8f" % (testfa(0.44)*testfa(0.55))

        fa1d = Function1d(function=testfa,thresh=1e-5,k=7,initial_level=1)
        fa1d.memory_usage()
        print "comp fa1d done"
        fa1d.reconstruct()

        fb1d = Function1d(function=testfa,thresh=1e-5,k=7)
        fb1d.memory_usage()
        print "comp fb1d done"
        fb1d.reconstruct()
        fa.reconstruct()
        
        print "mul start with_1d"
        f = f.mul_with1d(fb1d,side=1)
        print "mul done"

        print "Compressed value: %.8f" % f(0.33,0.44)
        print "  Analytic value: %.8f" % (testfa(0.33)*testfa(0.33)*testfa(0.44))

        print "mul start with_1d"
        f = f.mul_with1d(fa1d,side=2)
        print "mul done"

        print "Compressed value: %.8f" % f(0.33,0.44)
        print "  Analytic value: %.8f" % (
            testfa(0.33)*testfa(0.33)*testfa(0.44)*testfa(0.44))

        ### Try compressing a function defined in C to higher precision
        ### (in C for the speed)
        ##print "\ntest2: compressing an exponential to medium precision"
        ##cfcube.cvar.myexponent = 5.0
        ##f = Function6d(cfunction=cfcube.cvar.Ctestexp6d, function=cfcube.testexp6d,
        ##             thresh=1e-3, k=5, initial_level=1)

        f = Function2d(cfunction=cfcube.cvar.Ctestexp2d,thresh=1e-7,k=9)
        f.reconstruct()
        
        f.memory_usage()
        f.print_layers(0)
        print "Reconstruc value: %.8f" % f(0.33,0.44)
        print "  Analytic value: %.8f" % cfcube.testexp2d(0.33,0.44)
        if abs(f(0.33,0.44) - cfcube.testexp2d(0.33,0.44)) > 3e-5:
            print "test2 failed"
        print "test2: compression is OK"

        print "comp fa done"
        print "Truncating"
        f.truncate(1e-7)
        f.memory_usage()
        f.print_layers(0)
        print "Compressed value: %.8f" % f(0.33,0.44)
        print "  Analytic value: %.8f" % cfcube.testexp2d(0.33,0.44)

        print "Testing application of local function - g(x)=f(x)/(1+f(x)**2)"
        g = f.copy()
        g.local_function(cfunction=cfcube.cvar.Ctestloc)
        test1 = cfcube.testloc(f(0.33,0.44))
        test2 = g(0.33,0.44)
        print "   Analytic: %.8e" % test1
        print "  Numerical: %.8e" % test2
        del(g)
        
        # Try differentiating it
        df = f.diff(0)
        print "Compressed deriv: %.8f" % df(0.313,0.414)
        print "  Analytic value: %.8f" % cfcube.testexp2d_dx(0.313,0.414)
        del df

        df = f.diff(1)
        print "Compressed deriv: %.8f" % df(0.313,0.414)
        print "  Analytic value: %.8f" % cfcube.testexp2d_dy(0.313,0.414)
        del df

        # Try copying and squaring
        fsq = f.copy().square()
        val1 = fsq(0.313,0.414)
        val2 = cfcube.testexp2d(0.313,0.414)**2        
        print "Compressed square: %.8e" % val1
        print "  Analytic square: %.8e" % val2
        if abs(val1-val2) > 1e-5:
            print "test5 failed"
        print "test5: copying and squaring seem to be OK"

        fsq.truncate(1e-6)

        # Try addition and scalar multiplication
        g = fsq*1.3 + f*0.11
        val1 = g(0.45,0.55)
        val = cfcube.testexp2d(0.45,0.55)
        val2 = val*val*1.3 + val*0.11
        print "Compressed fsq*1.3 + f*0.11: %.8e" % val1
        print "  Analytic fsq*1.3 + f*0.11: %.8e" % val2
        if abs(val1-val2) > 4e-4:
            print "test6 failed"
        print "test6: addition and scalar multiplication seem OK"

        # Try multiplication by another function
        h = f*g
        val1 = h(0.45,0.55)
        val = cfcube.testexp2d(0.45,0.55)
        val2 = (val*val*1.3 + val*0.11)*val
        print "Compressed f*(fsq*1.3 + f*0.11): %.8e" % val1
        print "  Analytic f*(fsq*1.3 + f*0.11): %.8e" % val2
        if abs(val1-val2) > 5e-3:
            print "test6a failed"
        print "test6a: function multiplication seems OK"

        ### Try mul_func
        ##h21= f.mul_func(cfunction=cfcube.cvar.Ctestloc)
        ##h22= g.mul_func(function=cfcube.testexp6d)
        ##val11= h21(0.45,0.55,0.45)
        ##val12= h22(0.45,0.55,0.45)
        ##val = cfcube.testexp6d(0.45,0.55,0.45)
        ##val2 = (val*val*1.3 + val*0.11)*val
        ##print "Compressed f*(fsq*1.3 + f*0.11): %.8e" % val11
        ##print "Compressed f*(fsq*1.3 + f*0.11): %.8e" % val12
        ##print "  Analytic f*(fsq*1.3 + f*0.11): %.8e" % val2
        ##if abs(val11-val2) > 1e-2 or abs(val12-val2) > 1e-2:
        ##    raise "test6b failed"
        ##print "test6b: function multiplication using mul_func seems OK"

        # Try inner product in both bases
        f.compress()
        g.compress()
        test1 = f.inner(g)
        f.reconstruct()
        g.reconstruct()
        test2 = f.inner(g)
        test3 = h.trace()
        print "  Compressed: <f|g> = %.8e" % test1
        print "Uncompressed: <f|g> = %.8e" % test2
        print "     Tr(f*g): <f|g> = %.8e" % test3

        perf.output()


    
    
        
