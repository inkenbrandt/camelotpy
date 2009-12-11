#Copyright (C) 2009  Karl Bandilla
#
#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#at your option) any later version.
#
#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with this program.  If not, see http://www.gnu.org/licenses/.

import numpy as num
import scipy.special as sci
import math as math
import cmath as cmath
import HoogPy as hog

class SingleLayer:
    '''Implementation of solutions of head changes in aquifers
       due to injection based on
       Q. Zhou, J.T. Birkholzer, C.-F. Tsang
           (2009) A Semi-Analytical Solution for Large-Scale Perturbation
           and Leakage in a Laterally Bounded Aquifer-Aquitard System. Transport
           in Porous Media, 78(1), 127-148
       A.F. Moench (1985) Transient Flow to a Large-Diameter well in an aquifer
           with storative semi-confining layers. Water Resources Research,
           21(8), 1121-1131
       M.S. Hantush, C.E. Jacob (1955) Nonstedy radial flow in an infinite
           leaky aquifer. Trans. Am. Geophys. Union, 36th anual meeting,
           Pt.1, 95-100
       C.V. Theis (1935) The relation between lowering of piezometric surface
           and the rate and duration of discharge of a well using ground
           water storage. Trans. Am. Geophys. Union, 16th anual meeting,
           Pt.2, 519-524'''
    
    def __init__(self):
        #Set variables for Hoog transform
        self.N = 400
        self.work = num.zeros(self.N+1,dtype=num.complex128)
        self.d = num.zeros(self.N+1,dtype=num.complex128)
        self.meth = 3
        self.merror = 1e-10
                 
    def SetRequiredProperties(self, K, H, S, Q, Rwell, Rboundary = -1,\
                              LLboundarycond = 'constant', \
                              LLstorage = True, InjectionEnd = -1):
        '''Set the required minimum parameters to find pressure
           disturbance.
           K          list of hydraulic conductivities [L/T]
                      1 list entry for no leaky layers
                      2 list entries for aquifer + 1 leaky layer
                      3 list entries for aquifer + 2 leaky layers
                      aquifer value is the first one (S[0])
                      aquitard values must be in the same order as for
                      S and H
           H          list of thicknesses [L]
                      1 list entry for no leaky layers
                      2 list entries for aquifer + 1 leaky layer
                      3 list entries for aquifer + 2 leaky layers
                      aquifer value is the first one (S[0])
                      aquitard values must be in the same order as for
                      K and S
           S          list of storativites [-](volume of water that a
                      unit decline of head releases from storage in a
                      vertical prism of the aquifer of unit cross section)
                      1 list entry for no leaky layers
                      2 list entries for aquifer + 1 leaky layer
                      3 list entries for aquifer + 2 leaky layers
                      aquifer value is the first one (S[0])
                      aquitard values must be in the same order as for
                      K and H
           Q          injection rate [L^3/T] positve for injection
           Rwell      injection well radius [L]. Set to -1 for infinitesimal
                      well radius
           Rboundary  distance from injection well to lateral boundary [L]
                      default is set to -1, meaning no boundary (infinite
                      aquifer)
           LLboundarycond   boundary condition at the top of the upper
                      leaky layer and the bottom of the lower leaky layer
                      either 'constant' for constant head or 'noflow' for
                      no-flow boundary condition. 'constant' is set as default
           LLstorage  include leaky layer storage. either True or False.
                      True is set as default.
           InjectionEnd   time at which injection stopps [T]. Set to -1
                      for no injection end. Default set to -1'''
        #set conductivity list
        self.K = K
        #set thickness list
        self.H = H
        #set storativity list
        self.S = S
        #set injection rate
        self.Q = Q
        #set well radius
        self.rw = Rwell
        #set boundary radius
        self.rb = Rboundary
        #set leaky layer boundary condition
        self.LLbc = LLboundarycond
        #set if leaky layer storage is included
        self.LLstorage = LLstorage
        #set injection time
        self.InjectionEnd = InjectionEnd

    def SetProperty(self, proptype, value):
        '''Sets values for individual model properties
           proptype    'K': list of hydraulic conductivites
                       'H': list of layer thicknesses
                       'S': list of storativities
                       'Q': volumetric injection rate
                       'rw': well radius
                       'rb': boundary radius
                       'LLBC': leaky layer boundary condition
                       'LLstorage': is leaky layer storage included?
                       'InjectionEnd': end of injection time
           value       value or list of values to be assigned.'''
        #hydraulic conductivity
        if proptype == 'K':
            self.K = value
        #layer thickness
        if proptype == 'H':
            self.H = value
        #storativity
        if proptype == 'S':
            self.S = value
        #volumetric injection rate
        if proptype == 'Q':
            self.Q = value
        #well radius
        if proptype == 'rw':
            self.rw = value
        #boundary radius
        if proptype == 'rb':
            self.rb = value
        #leaky layer boundary condition
        if proptype == 'LLBC':
            self.LLbc = value
        #if leaky layer storage is included
        if proptype == 'LLstorage':
            self.LLstorage = value
        #end of injection time
        if proptype == 'InjectionEnd':
            self.InjectionEnd = value
            
    def SetDimensionlessParameters(self):
        '''Sets several dimensionless parameters needed for computation.
           SetRequiredParameters needs to be run before this function.'''
        #compute conversion factor for dimensionless time
        self.dimenTime = self.K[0] / self.S[0] / self.H[0]
        #set head factor
        self.HeadFactor = self.Q / (4.0 * math.pi * self.K[0] * self.H[0])
        #dimensionless well radius
        self.rDw = self.rw / self.H[0]
        #dimensionless boundary radius
        if self.rb > 0:
            self.rDb = self.rb / self.H[0]
        #compute leakage factors (lambda, equation (4f) from Zhou et al 2009)
        #storage factor (sigma, equation (5e) from Zhou et al 2009),
        #and m (equation (6d) from Zhou et al 2009) for the top aquitard
        #for the cases with leaky layers
        if len(self.S) > 1:
            #aquifer + 1 leaky layer
            #top aquitard
            self.lambdat = math.sqrt((self.K[1] / self.H[1]) \
                              / (self.K[0] / self.H[0]))
            if self.LLstorage:
                self.sigmat = self.S[1] / self.S[0]
            else:
                self.sigmat = 0.0
            self.mt = math.sqrt(self.sigmat) / self.lambdat
        if len(self.S) > 2:
            #aquifer + 2 leaky layers
            #bottom aquitard
            self.lambdab = math.sqrt((self.K[2] / self.H[2]) \
                              / (self.K[0] / self.H[0]))
            if self.LLstorage:
                self.sigmab = self.S[2] / self.S[0]
            else:
                self.sigmat = 0.0
            self.mb = math.sqrt(self.sigmab) / self.lambdab

    def GetHeadForRadiusTime(self, radius, time):
        '''Compute head using the solution by Zhou et al 2009
           radius   distance to injection well [L]
           time     time since injection started [T]'''
        #set head to zero if beyond boundary
        if self.rb > 0:
            if radius > self.rb:
                head = 0.0
                return head
        tD = float(time) * self.dimenTime
        self.rD = radius / self.H[0]
        bigT = tD * 3.0
        gamma = -math.log(self.merror) / (2.0 * bigT)
        init = 1
        head = hog.HoogTransform(tD, gamma, bigT, self.N, \
                                 self.LaplaceHead, self.meth, init, \
                                 self.d, self.work) * self.HeadFactor
        if self.InjectionEnd != -1 and time > self.InjectionEnd:
            tD = float(time - self.InjectionEnd) * self.dimenTime
            bigT = tD * 3.0
            gamma = -math.log(self.merror) / (2.0 * bigT)
            head = head - hog.HoogTransform(tD, gamma, bigT, self.N, \
                                 self.LaplaceHead, self.meth, init, \
                                 self.d, self.work) * self.HeadFactor
        return head

    def GetRadiusForHeadTime(self, head, time):
        '''Computes the contour radius from the injection well for a given
           head and time
           head   head for which contour is sought [L]
           time   time at which contour is sought [T]'''
        #set time variables
        tD = float(time) * self.dimenTime
        bigT = tD * 3.0
        gamma = -math.log(self.merror) / (2.0 * bigT)
        #start at the well radius, because highest head is here
        r1 = self.rw
        self.rD = r1 / self.H[0]
        init = 1
        res1 = hog.HoogTransform(tD, gamma, bigT, self.N, \
                                 self.LaplaceHead, self.meth, init, \
                                 self.d, self.work) * self.HeadFactor
        #find head at boundary (if exists), because there it will
        #be lowest
        if self.rb > 0:
            r2 = self.rb
            self.rD = r2 / self.H[0]
            init = 1
            res2 = hog.HoogTransform(tD, gamma, bigT, self.N, \
                                 self.LaplaceHead, self.meth, init, \
                                 self.d, self.work) * self.HeadFactor
        #if the head at the well radius is less than the specified
        #head, return 'no contour'
        if res1 < head:
            rcontour = 'no contour'
        #if the head at the boundary radius is greater than the specified
        #head, return 'no contour'
        elif self.rb > 0 and res2 > head:
            rcontour = 'no contour'
        else:
            #refine values are used to iteratively get closer to the
            #contour pressure
            refine = [0.3, 0.1, 0.05, 0.01, 0.01, 0.001]
            #multip values are used to bracket in the contour. so,
            #the first bracket is the well radius and 10 * the well radius.
            multip = [10, 50, 100, 500, 1000, 5000, 10000, 50000]
            #find the the bracket of the contour location
            for mul in multip:
                r2 = self.rw * mul
                #set r2 to rb if beyond boundary
                if self.rb > 0 and r2 > self.rb:
                    r2 = self.rb
                self.rD = r2 / self.H[0]
                init = 1
                res2 = hog.HoogTransform(tD, gamma, bigT, self.N, \
                                 self.LaplaceHead, self.meth, init, \
                                 self.d, self.work) * self.HeadFactor
                if res2 < head:
                    #if the pressure at r2 is less than the specified 
                    #countour the location bracket is found
                    break
                else:
                    #otherwise, the current result becomes the last result
                    #and the next location is tested
                    res1 = res2
                    r1 = r2
            #use linear interpolation to find contour location
            #log(0) is negative infinity, so set res2 to very small if 0.0
            if res2 == 0.0:
                res2 = 1.0e-12
            rcontour = self.Interpolate(r1, math.log(res1), r2, \
                                       math.log(res2), math.log(head))
            #the contour location is refined by finding the pressure
            #left and right from the the proposed contour and then
            #interpolating between the two values. the distance from
            #the proposed contour to left and right gets successively
            #smaller
            for dist in refine:
                #distance slightly smaller than proposed contour
                r1 = rcontour * (1.0 - dist)
                if r1 < self.rw:
                    r1 = self.rw
                self.rD = r1 / self.H[0]
                init = 1
                res1 = hog.HoogTransform(tD, gamma, bigT, self.N, \
                         self.LaplaceHead, self.meth, init, \
                         self.d, self.work) * self.HeadFactor
                #distance slightly larger than proposed contour
                r2 = rcontour * (1.0 + dist)
                if self.rb > 0 and r2 > self.rb:
                    r2 = self.rb
                self.rD = r2 / self.H[0]
                init = 1
                res2 = hog.HoogTransform(tD, gamma, bigT, self.N, \
                                 self.LaplaceHead, self.meth, init, \
                                 self.d, self.work) * self.HeadFactor
                #the new contour location is found by interpolation
                rcontour = self.Interpolate(r1, math.log(res1), r2, \
                                           math.log(res2), math.log(head))
        return rcontour

    def GetTimeForRadiusHead(self, radius, head):
        '''Compute time it takes to reach a certain head value at a given
           radius
           radius   distance from well at which head is sought [L]
           head     given head [L]'''
        self.rD = radius / self.H[0]

        #start at an arbitrary time
        tD = 1.0 * self.dimenTime
        old_res = -1.0
        
        #this captures any nan issues for large distances at early times
        while not old_res > 0.0:
            tD = 10.0 * tD
            bigT = tD * 3.0
            gamma = -math.log(self.merror) / (2.0 * bigT)
            init = 1
            old_res = hog.HoogTransform(tD, gamma, bigT, \
                            self.N, self.LaplaceHead, self.meth, \
                            init, self.d, self.work) * self.HeadFactor

        #if initial guess was too early increase time
        if old_res < head:
            #while the computed head is less than the given head increase the
            #elapsed time by a factor if 10 and recompute head
            while head > old_res:
                tD = 10.0 * tD
                if tD > 1e20:
                    return 'pressure not reached'
                bigT = tD * 3.0
                gamma = -math.log(self.merror) / (2.0 * bigT)
                init = 1
                old_res = hog.HoogTransform(tD, gamma, bigT, \
                            self.N, self.LaplaceHead, self.meth, \
                            init, self.d, self.work) * self.HeadFactor
            direction = 1.0
        #if initial guess was too late decrease time
        else:
            #while the computed head is greater than the given head
            #decrease the elapsed time by a factor if 10 and recompute
            #head
            while head < old_res:
                tD = tD / 10.0
                bigT = tD * 3.0
                gamma = -math.log(self.merror) / (2.0 * bigT)
                init = 1
                old_res = hog.HoogTransform(tD, gamma, bigT, \
                            self.N, self.LaplaceHead, self.meth, \
                            init, self.d, self.work) * self.HeadFactor
            direction = -1.0
        tstep = tD
        #once the computed head is greater than the given head, the time
        #is refined by going back-and-forth past the given head using ever
        #smaller time steps
        for i in range(0,5):
            #change the direction either from increasing time (direction=1)
            #to decreasing time (direction=-1) or vice-versa
            direction = direction * -1.0
            #decrease the time step by a factor of 10
            tstep = tstep / 10.0
            #for decreasing time the computed head starts higher than the
            #given time
            if direction == -1.0:
                higher = old_res
                lower = head
            #for increasing time the computed head starts lower than the
            #given time
            else:
                higher = head
                lower = old_res
            #as long as the higher head value remains higher change time by
            #one time step and recompute head
            while higher > lower:
                tD = tD + direction * tstep
                bigT = tD * 3.0
                gamma = -math.log(self.merror) / (2.0 * bigT)
                init = 1
                old_res = hog.HoogTransform(tD, gamma, bigT, \
                        self.N, self.LaplaceHead, self.meth, \
                        init, self.d, self.work) * self.HeadFactor
                #save recomputed head according to direction
                if direction == -1.0:
                    higher = old_res
                else:
                    lower = old_res
        #convert dimensionless time back to real time
        return tD / self.dimenTime
        
    def LaplaceHead(self, P):
        '''Laplace transform of dimensionless head (hD_bar) as given
           in its full form in equation (11) of Zhou et al 2009
           P    Laplace variable'''
        #f_bar as used in Zhou et al 2009 is the combination of the upper
        #and lower contributions
        #aquifer only, no leaky layer
        f_bar = 0.0
        #aquifer + 1 leaky layer
        if len(self.S) > 1:
            #from equation (6d) of Zhou et al 2009
            #mt is computed in the function FindZhouVariables
            fullmt = self.mt * cmath.sqrt(P)
            #the expression for f_bar in equation (8) of Zhou et al 2009
            #is simplified if leaky layer storage is neglected, or
            #the storage factor is very low
            if self.sigmat <= 2e-10:
                #simplified
                f_bar = self.lambdat**2
            else:
                #full expression
                if fullmt.real > 500.0:
                    f_bar = self.lambdat**2 * fullmt
                else:
                    #constant head leaky layer boundary condition
                    if self.LLbc == 'constant':
                        f_bar = self.lambdat**2 * fullmt \
                                / cmath.tanh(fullmt)
                    #no-flow leaky layer boundary condition
                    else:
                        f_bar = self.lambdat**2 * fullmt \
                                * cmath.tanh(fullmt)
        #aquifer + 2 leaky layers
        if len(self.S) > 2:
            #from equation (6d) of Zhou et al 2009
            #mb is computed in the function SetLBLVariables
            fullmb = self.mb * cmath.sqrt(P)
            #the expression for f_bar in equation (8) of Zhou et al 2009
            #is simplified if leaky layer storage is neglected, or
            #the storage factor is very low
            if self.sigmab <= 2e-10:
                #simplified
                f_bar += self.lambdab**2
            else:
                #full expression
                if fullmb.real > 500.0:
                    f_bar += self.lambdab**2 * fullmb
                else:
                    #constant head leaky layer boundary condition
                    if self.LLbc == 'constant':
                        f_bar += self.lambdab**2 * fullmb \
                                 / cmath.tanh(fullmb)
                    #no-flow leaky layer boundary condition
                    else:
                        f_bar = self.lambdab**2 * fullmb \
                                * cmath.tanh(fullmb)
        #from equation (10) of Zhou et al 2009
        x = cmath.sqrt(P + f_bar)
        #bounded aquifer
        if self.rb > 0:
            #to avoid numerical problems at early times
            temp = self.rDb * x
            #infinitesimal well radius
            if self.rw == -1:
                if temp.real > 500.0:
                    hD_bar = 2.0 / P  \
                       * sci.kv(0.0, self.rD * x)
                else: 
                    hD_bar = 2.0 / P  \
                       * (sci.kv(1.0, self.rDb * x) \
                          * sci.iv(0.0, self.rD * x) \
                       + sci.iv(1.0, self.rDb * x) \
                          * sci.kv(0.0, self.rD * x)) \
                                 / sci.iv(1.0, self.rDb * x)
            #finite well radius
            else:
                if temp.real > 500.0:
                    #to avoid numerical problems for VERY early times
                    #and close to the injection well
                    if sci.kv(0.0, x * self.rD) == 0j \
                               and sci.kv(1.0, self.rDw * x) == 0j:
                        hD_bar = 2.0 / (self.rDw * x * P)
                    else:
                        hD_bar = 2.0 * sci.kv(0.0, x * self.rD) \
                              / (self.rDw * x * P \
                                 * sci.kv(1.0, self.rDw * x))
                else:
                    hD_bar = 2.0 / (self.rDw * P * x) * \
                     (sci.kv(1.0, self.rDb * x) * sci.iv(0.0, self.rD * x) \
                    + sci.iv(1.0, self.rDb * x) * sci.kv(0.0, self.rD * x)) \
                   / (sci.iv(1.0, self.rDb * x) * sci.kv(1.0, self.rDw * x) \
                    - sci.kv(1.0, self.rDb * x) * sci.iv(1.0, self.rDw * x))
        #infinite aquifer
        else:
            #infinitesimal well radius
            if self.rw == -1:
                hD_bar = 2.0 / P * sci.kv(1.0, self.rDb * x)
            #finite well radius
            else:
                #to avoid numerical problems for VERY early times
                #and close to the injection well
                if sci.kv(0.0, x * self.rD) == 0j \
                           and sci.kv(1.0, self.rDw * x) == 0j:
                    hD_bar = 2.0 / (self.rDw * x * P)
                else:
                    hD_bar = 2.0 / (self.rDw * P * x) \
                         * sci.kv(0.0, self.rD * x) \
                         / sci.kv(1.0, self.rDw * x)
        return hD_bar

    def Interpolate(self, x1, y1, x2, y2, yknown):
        '''Linear interpolation between two points (x1,y1) and (x2,y2).
           yknown is the point between points 1 and 2 for which the
           interpolated value is sought.'''
        delx = x2 - x1
        dely = y1- y2
        return x2 - (yknown - y2) * delx / dely
