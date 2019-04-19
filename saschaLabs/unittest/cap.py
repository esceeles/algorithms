def capitalizer(text):
   """Capitalizes the first letter of each word in a string."""
   string = text.split()
   ans = ""
   for word in string:
      if word[0].isalpha():
         ans = ans + " " + word.capitalize()
      else:
         ans = ans + " " + word[0] + word[1:].capitalize()

   ans = ans.strip(" ")
   return ans
