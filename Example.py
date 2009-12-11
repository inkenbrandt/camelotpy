import SingleLayer_0_1_1 as sl
import Utilities as utl

# model properties
Ha = 500.0 #m, aquifer thickness
Ht = 200.0 #m, thickness of top leaky layer
Hb = 100.0 #m, thickness of bottom leaky layer
Q = 0.1 #m^3/s, injection rate
rw = 1.0 #m, well radius
rb = 60000.0 #m, boundary radius 
Ka = 1e-5 #m/s, aquifer hydraulic conductivity
Kt = 4e-12 #m2, hydraulic conductivity of top leaky layer
Kb = 4e-11 #m2, hydraulic conductivity of bottom leaky layer
Sa = 2e-3 #storativity of aquifer
St = 8e-4 #storativity of top leaky layer
Sb = 4e-4 #storativity of bottom leaky layer

#Specify when to generate results (time from injection start [s])
times = [1e2, 1e4, 1e6, 1e8, 1e9]
#Specify what pressure head [m] to generate radii for
contours = [0.1, 1.0, 10.0, 25.0]
#Specify distance from well [m]
radius = [10.0, 50.0, 100.0, 1000.0, 5000.0, 10000.0, 50000.0]

#initialize the model
model = sl.SingleLayer()
model.SetRequiredProperties([Ka, Kt, Kb], [Ha, Ht, Hb],  \
                            [Sa, St, Sb],Q, rw, rb,  \
                            'constant', True, -1)
model.SetDimensionlessParameters()

#Generate pressures for specified locations and times
print 'time [s], distance [m], head [m] in bounded domain'
for tim in times:
    for dis in radius:
        print tim, dis, model.GetHeadForRadiusTime(dis, tim)

#Generate pressure contours for specified pressures and times
print 'time [s], distance [m], head [m] in bounded domain'
for tim in times:
    for con in contours:
        print tim, model.GetRadiusForHeadTime(con, tim), con

#Generate times to reach pressure at specified distance
print 'time [s], distance [m], head [m] in bounded domain'
for rad in radius:
    for con in contours:
        print model.GetTimeForRadiusHead(rad, con), rad, con

#change model to infinite domain                             
model.SetProperty('rb', -1)
model.SetDimensionlessParameters()
#Generate pressures for specified locations and times
print 'time [s], distance [m], head [m] in an infinite domain'
for tim in times:
    for dis in radius:
        print tim, dis, model.GetHeadForRadiusTime(dis, tim)

#change change leaky layer boundary conditions                             
model.SetProperty('rb', rb)
model.SetProperty('LLBC', 'noflow')
model.SetDimensionlessParameters()
#Generate pressures for specified locations and times
print 'time [s], distance [m], head [m] in a bounded domain'
for tim in times:
    for dis in radius:
        print tim, dis, model.GetHeadForRadiusTime(dis, tim)

#end injection after t=1e7
model.SetProperty('LLBC', 'constant')
model.SetProperty('InjectionEnd', 1e7)
model.SetDimensionlessParameters()
#Generate pressures for specified locations and times
print 'time [s], distance [m], head [m] in an infinite domain'
for tim in times:
    for dis in radius:
        print tim, dis, model.GetHeadForRadiusTime(dis, tim)

