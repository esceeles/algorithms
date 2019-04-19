import itertools

file = open("uuid.txt")

#my hash function for part 1

print("Hashing uuids...")
hsh = {line.strip('\n'):sum([ord(l)**(i+1) for i,l in enumerate(line)]) for line in file}

print("collision rate: %{}".format(1-(len(set(hsh.values()))/len(hsh))))

#pt 2, using Python's built in hash function

class Hash_Map():
   def __init__(self):
      self.assoc_array = []
      for i in range (0, 1024):
         self.assoc_array.append(list())
      self.collisions = 0

   def Add(self, line):
      line = line.split(':')
      key = line[0].strip('"')
      value = line[1].strip('"')+":"+line[2].strip('"')+":"+line[3].strip('"')
      index = hash(key) % 1024
      if index in self.assoc_array:
         self.assoc_array[index].append(key, value)
         self.collisions += 1
      else:
         self.assoc_array[index] = (key, value)

   def In(self, key):
      for i in self.assoc_array:
         if key in i:
            return True
      return False

file.close()

file = open("key_value.txt")

a = Hash_Map()

#using file to fill in Hash Map
print("Building Hash Map...")
for line in file:
   a.Add(line.strip('\n'))

print(a.assoc_array)

print("collisions in Hash Map: ", a.collisions);
print("Is key \"Portland State University\" in Hash Map: ", a.In('Portland State University'))

file.close()
