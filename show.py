import numpy, os, tempfile,sys
import joblib 
def load(filename):
    return joblib.load(filename)["data"]

try:data = load(sys.argv[1])
except:data = load("block%s.dat" % sys.argv[1])
      
#if you want to change positions of the spheres along each segment, change these numbers
#e.g. [0,.1, .2 ...  .9] will draw 10 spheres, and this will look better
shifts = [0.,0.2,0.4,0.6,0.8]

#determining the 95 percentile distance between particles,  
meandist = numpy.percentile(numpy.sqrt(numpy.sum(numpy.diff(data,axis = 0)**2,axis = 1)),95)
#rescaling the data, so that bonds are of the order of 1. This is because rasmol spheres are of the fixed diameter. 
data /= meandist

#Screw you guys

#writing the rasmol script. Spacefill controls radius of the sphere. 
rascript = tempfile.NamedTemporaryFile()
rascript.write("""wireframe off 
color temperature
spacefill 100 
background white
""")
rascript.flush()


#creating the array, linearly chanhing from -225 to 225, to serve as an array of colors 
#(rasmol color space is -250 to 250, but it  still sets blue to the minimum color it found and red to the maximum). 
colors = numpy.array([int((j*450.)/(len(data)))-225 for j in xrange(len(data))])    

#creating spheres along the trajectory
#for speedup I just create a Nx4 array, where first three columns are coordinates, and fourth is the color      
newData = numpy.zeros((len(data) * len(shifts) - (len(shifts) - 1) ,4))  
for i in xrange(len(shifts)):            
    #filling in the array like 0,5,10,15; then 1,6,11,16; then 2,7,12,17, etc. 
    #this is just very fast
    newData[i:-1:len(shifts),:3] = data[:-1] * shifts[i] + data[1:] * ( 1 - shifts[i])            
    newData[i:-1:len(shifts),3] = colors[:-1]
newData[-1,:3] = data[-1]
newData[-1,3] = colors[-1]
            
towrite = tempfile.NamedTemporaryFile()
towrite.write("%d\n\n"%(len(newData)))  #number of atoms and a blank line after is a requirement of rasmol
    
for i in newData:                     
    towrite.write("CA\t%lf\t%lf\t%lf\t%d\n" % tuple(i)) 
towrite.flush()
#For windows you might need to change the place where your rasmol file is  
if os.name == "posix":  #if linux 
    os.system("rasmol -xyz %s -script %s" % (towrite.name, rascript.name))
else:     #if windows 
    os.system("C:/RasWin/raswin.exe -xyz %s -script %s" % (towrite.name, rascript.name))
    


#show3D(numpy.cumsum(numpy.random.randint(-1,2,(3,10000)),axis = 1))  #an example 

