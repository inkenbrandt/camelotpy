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

import math as mth
import scipy as spy
import scipy.special as sci

class HantushJacob:
    '''Implementation of solution of head change in aquifers due to
       injection based on the Hantush-Jacob solution: M.S. Hantush,
       C.E. Jacob (1955) Nonsteady radial flow in an infinite leaky aquifer.
       Transactions of the American Geophysical Union, 1, 95-100'''
    def __init__(self):
        pass

    def FindB(self, aquifer_transmissivity, aquitard_conductivities, \
              aquitard_thicknesses):
        '''Finds B parameter for Hantush-Jacob solution.
           aquifer_transmissivity   aquifer transmissivity [L^2/T] 
           aquitard_conductivities  list of aquitard hydraulic conductvities
                                    minimum 1 entry, maximum 2 entries [L/T]
                                    needs to have same number of entries as
                                    aquitard_thicknesses
           aquitard_thicknesses     list of aquitard thicknesses
                                    minimum 1 entry, maximum 2 entries [L]
                                    needs to have same number of entries as
                                    aquitard_conductivities'''
        #find length of aquitard lists
        Naquitards = len(aquitard_conductivities)
        #for one aquitard
        if Naquitards == 1:
            return mth.sqrt(aquifer_transmissivity \
                            * aquitard_thicknesses[0] \
                            / aquitard_conductivities[0])
        #for two aquitards
        elif Naquitards == 2:
            return mth.sqrt(aquifer_transmissivity \
                            * aquitard_thicknesses[0] \
                            * aquitard_thicknesses[1] \
                            /(aquitard_conductivities[0] \
                              * aquitard_thicknesses[1] \
                              + aquitard_conductivities[1] \
                              * aquitard_thicknesses[0]))
        
            
    def SetRequiredProperties(self, Q, S, T, B, InjectionEnd = -1, \
                              HeadOrPressure = 'H', FreshWDensity = 1000.0, \
                              Gravity = 9.81, WarnNoRadius=["s", \
                              "no contour"], WarnGTsteadystate=["s", \
                              "pressure not reached"]):
        '''Set the required minimum parameters to find head
           disturbance.
           Q       injection rate [L^3/T] positve for injection
           S       storativity of the aquifer (volume of water that a
                   unit decline of head releases from storage in a
                   vertical prism of the aquifer of unit cross section) [-]
           T       transmissivity of aquifer [L^2/T]
           B       square root of transmissivity divided by leakance
           InjectionEnd   time at which injection stopps [T]. Set to -1
                   for no injection end. Default set to -1
           HeadOrPressure   should the output be pressure ("P") or
                            piezometric ("H")
           FreshWDensity    fresh water density, is necessary if pressure 
                            and not head is to be computed [M/L^3]
           Gravity    gravitational accleration [L/T^2] is necessary if
                      pressure and not head is to be computed
           WarnNoRadius   value to return if no distance can be
                          found in GETRadiusForHeadTime. First value is
                          either "s" for a string or "i" for an integer.
                          Second value is returned.
           WarnGTsteadystate   value to return if no time can be
                          found in GetTimeForRadiusHead, becuase head is
                          bellow staeady state head. First value is
                          either "s" for a string or "i" for an integer.
                          Second value is returned.'''
        #set injection rate
        self.Q = Q
        #set storativity
        self.S = S
        #set transmissivity
        self.T = T
        #set B
        self.B = B
        #set injection time
        self.InjectionEnd = InjectionEnd
        #set head or pressure
        self.HorP = HeadOrPressure
        #set fresh water density
        self.FreshWDensity = FreshWDensity
        #set gravitational accleration
        self.Gravity = Gravity
        #set head factor
        self.HeadFactor = self.Q / (4.0 * mth.pi * self.T)
        if self.HorP == 'P':
            self.HeadFactor = self.HeadFactor * self.Gravity \
                              * self.FreshWDensity
        #set u factor
        self.uFactor = self.S / (4.0 * self.T)
        #set warning messages
        if WarnNoRadius[0] == "s":
            self.WarnNoRadius = str(WarnNoRadius[1])
        elif WarnNoRadius[0] == "i":
            self.WarnNoRadius = float(WarnNoRadius[1])
        if WarnGTsteadystate[0] == "s":
            self.WarnGTsteadystate = str(WarnGTsteadystate[1])
        elif WarnGTsteadystate[0] == "i":
            self.WarnGTsteadystate = float(WarnGTsteadystate[1])

    def SetProperty(self, proptype, value):
        '''Sets values for individual model properties
           proptype    'T': aquifer transmissivity
                       'S': aquifer storativity
                       'Q': volumetric injection rate
                       'B': transmissivity-leakance ratio
                       'InjectionEnd': end of injection time
                       'HeadOrPressure': output either head or pressure
                       'FreshWDensity': fresh water density
                       'Gravity': gravity
                       'WarnNoRadius': return value if no radius is found
                       'WarnGTsteadystate': return value if pressure is
                                            above steady state pressure
           value       value or list of values to be assigned.'''
        #hydraulic conductivity
        if proptype == 'T':
            self.T = value
        #storativity
        if proptype == 'S':
            self.S = value
        #volumetric injection rate
        if proptype == 'Q':
            self.Q = value
        #transmissivity-leakance ratio
        if proptype == 'B':
            self.B = value
        #end of injection time
        if proptype == 'InjectionEnd':
            self.InjectionEnd = value
        #output head or pressure
        if proptype == 'HeadOrPressure':
            self.HorP = value
        #brine density
        if proptype == 'FreshWDensity':
            self.FreshWDensity = value
        #gravitational accleration
        if proptype == 'Gravity':
            self.Gravity = value
        #no radius warning
        if proptype == 'WarnNoRadius':
            if value[0] == "s":
                self.WarnNoRadius = str(value[1])
            elif value[0] == "i":
                self.WarnNoRadius = float(value[1])
        #pressure not reached warning
        if proptype == 'WarnGTsteadystate':
            if value[0] == "s":
                self.WarnGTsteadystate = str(value[1])
            elif value[0] == "i":
                self.WarnGTsteadystate = float(value[1])
        #set HJ factor
        self.HeadFactor = self.Q / (4.0 * mth.pi * self.T)
        if self.HorP == 'P':
            self.HeadFactor = self.HeadFactor * self.Gravity \
                              * self.FreshWDensity
        #set u factor
        self.uFactor = self.S / (4.0 * self.T)

    def EvalWellFunction(self, u, rB):
        '''Evaluates the Hantush-Jacob well function.
           u     dimensionless time
           rB    dimensionless radius'''
        #return zero for large u to avoid numerical issues
        if u > 14:
            return 0.0        
        #set number of terms used in summation. as u increases more
        #terms are needed
        if u <= 2:
            endmember = 10
        elif u <= 7:
            endmember = 20
        else:
            endmember = 30
        #compute temporary variable
        r4B = rB**2 / 4.0
        #compute summation term
        last_term = 0.0
        for i in range(1,endmember):
            for j in range(1,i+1):
                last_term += (-1)**(i+j) * spy.factorial(i - j + 1) / \
                             (spy.factorial(i+2))**2 * u**(i-j) * r4B**j
        #for rB values close to 0 one term is dropped to avoid numerical
        #problems
        if abs(sci.iv(0.0,rB) - 1.0) < 1e-15:
            WellFuncHJ = 2.0 * sci.kv(0.0,rB) \
                 - sci.iv(0.0,rB) * sci.expn(1,r4B / u) \
                 + mth.exp(-1.0 * r4B / u) * (0.5772156649015328606 \
                 + mth.log(u) + sci.expn(1,u) - u - u**2 * last_term)
        else:
        #full expression
            WellFuncHJ = 2.0 * sci.kv(0.0,rB) \
                 - sci.iv(0.0,rB) * sci.expn(1,r4B / u) \
                 + mth.exp(-1.0 * r4B / u) * (0.5772156649015328606 \
                 + mth.log(u) + sci.expn(1,u) - u + u * (sci.iv(0.0,rB) - 1.0) \
                    / r4B - u**2 * last_term)
        return WellFuncHJ
        
            
    def ListOrSingle(self, inValues):
        '''Determines if the inValues are a list or a single value.
           Returns inValues in the form of a list and a boolean that
           is True if tinValues are a list and False if innValues is
           a single value.'''
        if isinstance(inValues, list):
            outList = inValues
            ListOut = True
        else:
            outList = [inValues]
            ListOut = False
        return outList, ListOut

    def GetHeadForRadiusTime(self, radius, time):
        '''Compute head using the Hantush-Jacob solution
           radius   distance to injection well [L]
           time     time since injection started [T]'''
        #check if inputs are scalars or lists
        RadiusList, RadIsList = self.ListOrSingle(radius)
        TimeList, TimIsList = self.ListOrSingle(time)
        if not RadIsList and not TimIsList:
            ListOut = False
        else:
            ListOut = True
        #initialize output list
        Heads = []
        for tim in TimeList:
            for rad in RadiusList:
                #set increase to zero for time = 0.0
                if tim == 0.0:
                    Heads.append([tim, rad, 0.0])
                else:
                    #find u
                    u = self.uFactor * rad**2 / tim
                    #set increase to zero for large u
                    if u > 14:
                        Heads.append([tim, rad, 0.0])
                    else:
                        head = self.HeadFactor \
                               * self.EvalWellFunction(u, rad / self.B)
                        #for times after injection has ceased
                        if self.InjectionEnd != -1 \
                                        and tim > self.InjectionEnd:
                            #find u
                            u = self.uFactor * rad**2 \
                                        / (tim - self.InjectionEnd)
                            #set increase to zero for large u
                            if u <= 14:            
                                head = head - self.HeadFactor \
                                        * self.EvalWellFunction(u, \
                                                        rad / self.B)
                            if head < 0.0:
                                head = 0.0
                        Heads.append([tim, rad, head])
        if ListOut:
            return Heads
        else:
            return Heads[0][2]


    def GetRadiusForHeadTime(self, head, time):
        '''Computes the contour radius from the injection well for a given
           head and time from Theis solution
           head   head for which contour is sought [L]
           time   time at which contour is sought [T]'''
        #check if inputs are scalars or lists
        HeadList, HedIsList = self.ListOrSingle(head)
        TimeList, TimIsList = self.ListOrSingle(time)
        if not HedIsList and not TimIsList:
            ListOut = False
        else:
            ListOut = True
        #initialize output list
        Distances = []
        for tim in TimeList:
            for hed in HeadList:
                #start at the well radius, because highest head is here
                smallest_r = 0.0001
                r1 = smallest_r
                res1 = self.GetHeadForRadiusTime(r1, tim)
                #if the head at the well radius is less than the specified
                #head, return 'no contour'
                if res1 < hed:
                    Distances.append([tim, self.WarnNoRadius, hed])
                else:
                    #refine values are used to iteratively get closer to the
                    #contour pressure
                    refine = [0.3, 0.1, 0.05, 0.01, 0.01, 0.001]
                    #multip values are used to bracket in the contour. so,
                    #the first bracket is smallets_r and 10 * smallest_r.
                    multip = [10, 50, 100, 500, 1000, 5000, 10000, \
                              50000, 100000, 500000, 1000000, 5000000, \
                              10000000, 50000000, 100000000, 500000000, \
                              1000000000, 5000000000, 10000000000, \
                              50000000000, 100000000000, 500000000000]
                    #find the the bracket of the contour location
                    for mul in multip:
                        r2 = smallest_r * mul
                        res2 = self.GetHeadForRadiusTime(r2, tim)
                        if res2 < hed:
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
                    rcontour = self.Interpolate(r1, mth.log(res1), r2, \
                                                mth.log(res2), \
                                                mth.log(hed))
                    #the contour location is refined by finding the pressure
                    #left and right from the the proposed contour and then
                    #interpolating between the two values. the distance from
                    #the proposed contour to left and right gets successively
                    #smaller
                    for dist in refine:
                        #distance slightly smaller than proposed contour
                        r1 = rcontour * (1.0 - dist)
                        if r1 <= smallest_r:
                            r1 = smallest_r
                        res1 = self.GetHeadForRadiusTime(r1, tim)
                        #distance slightly larger than proposed contour
                        r2 = rcontour * (1.0 + dist)
                        res2 = self.GetHeadForRadiusTime(r2, tim)
                        #the new contour location is found by interpolation
                        rcontour = self.Interpolate(r1, mth.log(res1), \
                                                    r2, mth.log(res2), \
                                                    mth.log(hed))
                    Distances.append([tim, rcontour, hed])
        if ListOut:
            return Distances
        else:
            return Distances[0][1]

    def GetTimeForRadiusHead(self, radius, head):
        '''Compute time it takes to reach a certain head value at a given
           radius
           radius   distance from well at which head is sought [L]
           head     given head [L]'''
        #check if inputs are scalars or lists
        RadiusList, RadIsList = self.ListOrSingle(radius)
        HeadList, HedIsList = self.ListOrSingle(head)
        if not RadIsList and not HedIsList:
            ListOut = False
        else:
            ListOut = True
        #initialize output list
        Times = []
        for rad in RadiusList:
            for hed in HeadList:
                 #start at an arbitrary time
                time = 1.0
                old_res = -1.0
                #this captures any nan issues for large distances at early times
                while not old_res > 0.0:
                    time = time * 10.0
                    old_res = self.GetHeadForRadiusTime(rad, time)
       
                #outlet if head is not reached
                finished = False
                #if initial guess was too early increase time
                if old_res < hed:
                    #while the computed head is less than the given head increase the
                    #elapsed time by a factor if 10 and recompute head
                    while hed > old_res:
                        time = time * 10.0
                        if time > 3.15e10:
                            Times.append([self.WarnGTsteadystate, rad, \
                                          hed])
                            finished = True
                            break
                        old_res = self.GetHeadForRadiusTime(rad, time)
                    direction = 1.0
                #if initial guess was too late decrease time
                else:
                    #while the computed head is greater than the given head
                    #decrease the elapsed time by a factor if 10 and recompute
                    #head
                    while hed < old_res:
                        time = time /10.0
                        old_res = self.GetHeadForRadiusTime(rad, time)
                    direction = -1.0
                timestep = time
                if not finished:
                    #once the computed head is greater than the given head, the time
                    #is refined by going back-and-forth past the given head using ever
                    #smaller time steps
                    for i in range(0,5):
                        #change the direction either from increasing time (direction=1)
                        #to decreasing time (direction=-1) or vice-versa
                        direction = direction * -1.0
                        #decrease the time step by a factor of 10
                        timestep = timestep / 10.0
                        #for decreasing time the computed head starts higher than the
                        #given time
                        if direction == -1.0:
                            higher = old_res
                            lower = hed
                        #for increasing time the computed head starts lower than the
                        #given time
                        else:
                            higher = hed
                            lower = old_res
                        #as long as the higher head value remains higher change time by
                        #one time step and recompute head
                        while higher > lower:
                            time = time + direction * timestep
                            old_res = self.GetHeadForRadiusTime(rad, time)
                            #save recomputed head according to direction
                            if direction == -1.0:
                                higher = old_res
                            else:
                                lower = old_res
                    #convert dimensionless time back to real time
                    Times.append([time, rad, hed])
        if ListOut:
            return Times
        else:
            return Times[0][0]

    def Interpolate(self, x1, y1, x2, y2, yknown):
        '''Linear interpolation between two points (x1,y1) and (x2,y2).
           yknown is the point between points 1 and 2 for which the
           interpolated value is sought.'''
        delx = x2 - x1
        dely = y1- y2
        return x2 - (yknown - y2) * delx / dely
