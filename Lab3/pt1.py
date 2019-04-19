file = open("data5.txt", "r") #opens file in read mode
textCount = int(file.readline())   #gets first line for num of chars in text

file.readline()

text = file.read(textCount) #gets text

file.readline() #throws out newline 
file.readline()

patternCount = int(file.readline())     #gets number of chars in pattern

file.readline() 

pattern = file.read(patternCount)       #gets pattern

print("\nPattern:\n", pattern)
index = text.find(pattern)              #finds index of pattern in string

print("\nText:\n", text)

b = ""                  #creates string with proper formatting for carot placement
for i in range (0 , index):
    b += " "

print(b, "^")

print("\n")

file.close()   
