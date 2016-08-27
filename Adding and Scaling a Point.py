#/Users/Apollo/miniconda/bin/python
""" This problem will add two points together, or Scale a point by a scalar value """



x1 = int(raw_input('Enter your Point 1 x-coordinate :'))
y1 = int(raw_input('Enter your Point 1 y-coordinate :'))
point_1 = (x1,y1)
print "Point 1 is : %s" % (point_1,)
x2 = int(raw_input('Enter your Point 2 x-coordinate :'))
y2 = int(raw_input('Enter your Point 2 y-coordinate :'))
point_2 = (x2,y2)
print "Point 2 is : %s" % (point_2,)

class Point(object):
    def __init__(self, x, y):
        self.x = x
        self.y = y
        
        
P1 = Point(x1,y1)
P2 = Point(x2,y2)

P3 = (P1.x + P2.x, P1.y + P2.y)
print '---------------------------'
print '---------------------------'
print "Point 1 + Point 2 =  %s" % (P3,)

question = raw_input('Would you like to scale your point by a scalar i.e u*(x,y) ? [Y/N] : ')

if question == 'Y':
    u = int(raw_input('Enter your scalar value, u :'))
    question_2 = raw_input('Which Point would you like to scale ? [P1,P2,P3] :')
    if question_2 == 'P1':
        P1_new = (u*P1.x, u*P1.y)
        print P1_new
    elif question_2 == 'P2':
        P2_new = (u*P2.x, u*P2.y)
        print P2_new      
    elif question_2 == 'P3':
        P3_new = (u*(P1.x + P2.x), u*(P1.y + P2.y))
        print P3_new
    else:
        pass
