import numpy as np
import matplotlib.pyplot as plt

class Point(object) :
    
    def __init__(self, x, y) :
        self.x, self.y = x, y

    def __add__(self,other):
        x = self.x + other.x
        y = self.y + other.y
        return Point(x,y)
    
    def __mul__(self,scale):
        x = self.x * scale
        y = self.y * scale
        return Point(x,y)
    
    def __rmul__(self,other):
        return Point(self.x*other,self.y*other)
    
    def __str__(self) :
        return "Point(%.6F, %.6f) " % (self.x, self.y)

class Ray(object) :
    
    def __init__(self, origin, direction) :
        self.origin = origin
        # ensure the direction is normalized to unity, i.e., cos^2 + sin^2 = 1
        norm = np.sqrt(direction.x**2 + direction.y**2)
        self.direction = Point(direction.x/norm, direction.y/norm)
            
    def __str__(self) :
        return "Ray: r_0(%10.6f, %10.6f), d(%.6f %.6f) " % \
               (self.origin.x, self.origin.y, self.direction.x, self.direction.y)
   
class Node(object) :
    
    def contains(self, p) :
        """Does the node contain the point?"""
        raise NotImplementedError

    def intersections(self, r) :
        """Where does the node intersect the ray?"""
        raise NotImplementedError

class Primitive(Node) :

    def __init__(self, surface, sense) :
        self.surface, self.sense = surface, sense

    def contains(self, p) :
        return (self.surface.f(p) < 0) == self.sense

    def intersections(self, r) :
        return self.surface.intersections(r)

class Operator(Node) :
    
    def __init__(self, L, R) :
        self.L, self.R = L, R
        
    def contains(self, p) :
        raise NotImplementedError

    def intersections(self, r) :
        # get intersections with left and right nodes
        pointsL = self.L.intersections(r)
        pointsR = self.R.intersections(r)
        # return the concatenated result
        return pointsL + pointsR
             
class Union(Operator):
    
    def contains(self,p):
        
        if self.L.contains(p) or self.R.contains(p):
            return True
        else:
            return False
   
class Intersection(Operator):
        
    def contains(self,p):
        if self.L.contains(p) and self.R.contains(p):
            return True
        else:
            return False
            
class Surface(object) :
    
    def f(self, p) :
        raise NotImplementedError
        
    def intersections(self, r) :
        raise NotImplementedError
             
class QuadraticSurface(Surface) :
    
    def __init__(self, A=0.0, B=0.0, C=0.0, D=0.0, E=0.0, F=0.0) :
        super(QuadraticSurface,self).__init__()
        self.A = A
        self.B = B
        self.C = C
        self.D = D
        self.E = E
        self.F = F
    
    def intersections(self, r) :
        self.r = r
        M = np.array([[2*self.A, self.C, self.D],[self.C, 2*self.B, self.E],[self.D, self.E, 2*self.F]])
        r0 = np.array([r.origin.x,r.origin.y,1.0])
        d = np.array([r.direction.x,r.direction.y,0.0])
        a = np.dot(np.dot(np.transpose(d),M),d)
        b = 2*np.dot(np.dot(np.transpose(r0),M),d)
        c = np.dot(np.dot(np.transpose(r0),M),r0)
        undersqrt = (b**2) - 4*a*c
        
        if a == 0.0:
            t = -c/b
            p = r.origin + t*r.direction
            return [p]
        elif undersqrt > 0.0:
            t1 = (-b + np.sqrt(undersqrt))/(2*a)
            t2 = (-b - np.sqrt(undersqrt))/(2*a)
            p1 = r.origin + t1*r.direction
            p2 = r.origin + t2*r.direction
            return [p1,p2]
        elif undersqrt == 0.0:
            t = -b/2*a
            p = r.origin + t*r.direction
            return [p]
        elif undersqrt < 0.0:
            return []
        
    def f(self, p) :
        x = p.x
        y = p.y
        return self.A*x**2 + self.B*y**2 + self.C*x*y + self.D*x + self.E*y + self.F
               
class Region(object) :
    
    def __init__(self) :
        self.node = None
    # We use this method to put surfaces and nodes together.
    # This method is passed node, surface, and it's operation 'Union' or 'Intersection'
    # Once implemented, we can use the region.append() to append new surfaces
    # Enabling us to graph various surfaces we have defined already such as Circles, Planes etc.
    # As well as their Intersections and/or Unions
    def append(self, node=None, surface=None, operation="U", sense=False) :
        # This assertion is essentially an internal self check method
        # This will ensure when the code runs, it (hopefully) runs without errors        
        assert((node and not surface) or (surface and not node))
        # Determines if the instance is a surface
        # If the instance is indeed a surface, then assigns node to the class Primitive
        if isinstance(surface, Surface) :
            node = Primitive(surface, sense)
            
        if self.node is None :
            self.node = node
           
        else :
            O = Union if operation == "U" else Intersection
            self.node = O(self.node, node)
        
    # Orders points intersected by ray from first encountered to last
    def intersections(self, r) :
        ints = self.node.intersections(r)
        ints.sort(key = lambda p: p.x)
        return ints
        
    def contains(self,p):
        return self.node.contains(p)
             
class Geometry(object) :
    
    noregion = -1    
    newregion = []
    def __init__(self,  xmin, xmax, ymin, ymax) :
        self.xmin, self.xmax, self.ymin, self.ymax = xmin, xmax, ymin, ymax
        self.regions = []
          
    def add_region(self, r) :
        self.regions.append(r)
        
    def find_region(self, p) :
        region = Geometry.noregion
        newregion = []
        
        if p.x > self.xmax or p.x < self.xmin or p.y > self.ymax or p.y < self.ymin:
            newregion.append(self.noregion) 
        else:
            
            for i in self.regions:
                region = i
                if region.contains(p):
                    newregion.append(self.regions.index(region))
            if len(newregion) == 0:
                newregion.append(self.noregion)
            
            return newregion
        
    def plot(self, nx, ny) :
        x = np.linspace(self.xmin,self.xmax,nx)
        y = np.linspace(self.ymin,self.ymax,ny)
        
        # Generate random colors for regions
        colors = np.random.rand(len(self.regions))
        
        x_i = []
        y_i = []
        color_i = []
        
        for i in x:
            for j in y:
                stop = 0
                testPoint = Point(i,j)
                Region_i = self.find_region(testPoint)
                
                if Region_i[0] == -1:
                    stop = 0
                elif len(Region_i) > 1:
                    color = 'g'
                    stop = 1
                else :
                    color = colors[Region_i[0]]
                    stop = 1
                    
                if stop == 1:
                    x_i.append(i)
                    y_i.append(j)
                    color_i.append(color)
                    
        plt.scatter(x_i, y_i, c=color_i, marker ='.', edgecolors='None')
        
class PlaneV(QuadraticSurface):
    def __init__(self, a):
        super(PlaneV,self).__init__(D = 1.0, F=-a)

class PlaneH(QuadraticSurface):
    def __init__(self, b):
        super(PlaneH,self).__init__(A=0.0,B=0.0,C=0.0,D=0.0,E=0.0,F=0.0)
        self.b = b       
        self.A = 0.0
        self.B = 0.0
        self.C = 0.0
        self.D = 0.0
        self.E = 1.0
        self.F = -self.b

class Plane(QuadraticSurface):
    def __init__(self,m,b):
        self.m = m
        self.b = b
        self.A = 0.0
        self.B = 0.0
        self.C = 0.0
        self.D = -1.0
        self.E = 1.0
        self.F = -self.b        

class Circle(QuadraticSurface):
    def __init__(self,r,a=0.0,b=0.0):
        super(Circle,self).__init__(A=0.0,B=0.0,C=0.0,D=0.0,E=0.0,F=0.0)
        self.r = r
        self.a = a
        self.b = b
        self.A = 1.0
        self.B = 1.0
        self.C = 0.0
        self.D = -2.0*self.a
        self.E = -2.0*self.b
        self.F = self.a**2 + self.b**2 - self.r**2
        
# Constructs two circle regions, and plots them
c0 = Circle(1)
c1 = Circle(2)
region_new = Region()
region_new.append(surface = c0, operation='I', sense=False)
region_new.append(surface = c1, operation='I', sense=True)
region_new2 = Region()
region_new2.append(surface = c0, operation='U', sense=True)

G = Geometry(-10,10,-10,10)
G.add_region(region_new)
G.add_region(region_new2)

G.plot(500,500)



