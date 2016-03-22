import HantushJacob_0_2 as hj

# model properties
Ha = 500.0 #m, aquifer thickness
Ht = 200.0 #m, thickness of top leaky layer
Hb = 100.0 #m, thickness of bottom leaky layer
Q = 0.1 #m^3/s, injection rate
Ka = 1e-5 #m/s, aquifer hydraulic conductivity
Kt = 4e-12 #m2, hydraulic conductivity of top leaky layer
Kb = 4e-11 #m2, hydraulic conductivity of bottom leaky layer
S = 2e-3 #storativity of aquifer
T = 5e-3 #m2/s, aquifer transmissivity
#Specify when to generate results (time from injection start [s])
times = [1e2, 1e4, 1e6, 1e9]
#Specify what pressure head [m] to generate radii for
contours = [0.1, 10.0, 50.0]
#Specify distance from well [m]
radius = [10.0, 100.0, 5000.0, 10000.0, 100000.0]

#initialize the model
model = hj.HantushJacob()

#compute B
B = model.FindB(T, [Kt, Kb], [Ht, Hb])
model.SetRequiredProperties(Q, S, T, B)

#Generate heads for specified locations and times
print 'Head for given time and distance'
print 'time [s], distance [m], head [m]'
for tim in times:
    for dis in radius:
        print tim, dis, model.GetHeadForRadiusTime(dis, tim)

#Generate head contours for specified pressures and times
model.SetProperty('WarnNoRadius', ['i', -2])
print 'Contour radius for given head and time'
print 'time [s], distance [m], head [m]'
for tim in times:
    for con in contours:
        print tim, model.GetRadiusForHeadTime(con, tim), con

#Generate times to reach head at specified distance
print 'Time to reach given head at given radius'
print 'time [s], distance [m], head [m]'
for rad in radius:
    for con in contours:
        print model.GetTimeForRadiusHead(rad, con), rad, con

#Generate pressures for specified locations and times with an injection
#shut-off after 1e7 seconds
model.SetProperty('InjectionEnd', 1e7)
model.SetProperty('HeadOrPressure', 'P')
model.SetProperty('FreshWDensity', 1000.0)
print 'Head for given time and distance with shut-off after 1e7 s'
print 'using list input and output'
print 'time [s], distance [m], pressure [Pa]'
print model.GetHeadForRadiusTime(radius, times)
