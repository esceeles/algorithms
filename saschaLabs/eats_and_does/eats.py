#eats and does assignment for Charlie and Ellie
class Bear():
    def __init__(self, name):
        self.name = name
    def __lt__ (self, other):
        return self.name < other.name
    def __le__ (self, other):
        return self.name <= other.name
    def __gt__ (self, other):
        return self.name > other.name
    def __ge__ (self, other):
        return self.name >= other.name
    def eats(self):
        return 'berries'

class Turtle():
    def __init__(self, name, readingLevel):
        self.name = name
        self.level = readingLevel
    def eats(self):
        return 'insects and small fish'
    def __lt__ (self, other):
        return self.level < other.level
    def __le__ (self, other):
        return self.level <= other.level
    def __gt__ (self, other):
        return self.level > other.level
    def __ge__ (self, other):
        return self.level >= other.level

class Sasquatch():
    def __init__(self, name, areaCode):
        self.name = name
        self.ac = areaCode
    def eats(self):
        return 'hazelnuts'
    def __lt__ (self, other):
        return self.ac < other.ac
    def __le__ (self, other):
        return self.ac <= other.ac
    def __gt__ (self, other):
        return self.ac > other.ac
    def __ge__ (self, other):
        return self.ac >= other.ac


Baloo = Bear("Baloo")
Tammy = Turtle("Tammy", 4)
Smithy = Sasquatch("Smith", 218)

print(f'The bear named Baloo eats {Baloo.eats()}')
print(f'The turtle named Tammy eats {Tammy.eats()}')
print(f'The sasquatch named Smith eats {Smithy.eats()}')

Smokey = Bear("Smokey")
Timothy = Turtle("Timothy", 12)
Harry = Sasquatch("Harry", 503)

print("Is the name Smokey < Baloo?", Smokey < Baloo)
print("Is Tammy's reading level < Timothy's?", Tammy < Timothy)
print("Is Smithy's area code >= Harry's?", Smithy >= Harry)
