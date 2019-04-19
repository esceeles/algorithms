import requests


def create_lang_abbr_map():
    """
    This function creates a dictionary that maps language abbreviations,
    "codes", used by detectlanguage.com to full language names.
    """
    r = requests.get("https://ws.detectlanguage.com/0.2/languages")
    return {x["code"]: x["name"] for x in r.json()}


lang_abbr_map = create_lang_abbr_map()


def format_lang_string(lang_str):
    """
    This function formats a string to be a parameter in a POST request
    to https://detectlanguage.com
    """
    lang_str = lang_str.replace(' ', '+')
    lang_str = 'q=' + lang_str
    return lang_str


def post_language_samples():
    lang_samples = ["Përshendetje Botë", "Русский компьютер", "Ahoj světe",
                    "ਕੰਪਿਊਟਰ ਵਿਗਿਆਨ", "مرحبا بالعالم", "Labas pasauli",
                    "Привет мир", "ਯੇਅ ਨਿਊ ਸ਼ੁਰੂਆਤ", "ہیلو دنیا",
                    "Salam dünya", "नमस्कार संसार", "Halló heimur",
                    "Беларускія Кампутарныя", "ਸਤਿ ਸ੍ਰੀ ਅਕਾਲ ਦੁਨਿਆ",
                    "Прывітанне Сусвет"]

    r = list()    
    for i in lang_samples:
      x = requests.post("https://ws.detectlanguage.com/0.2/detect", headers= {"Authorization" : "Bearer 80b4ffd09ed8d997433f36a621dc9d62"}, params=(format_lang_string(i)))
      r.append(x)
    m = create_lang_abbr_map() 
    l = list()
    for i in r:
      l.append(get_lang_abbr_from_resp(i))
    seen = {}
    most = 0
    for i in l:
      if l.count(i) > most:
         x = i
         most = l.count(i)
    print("most frequent language is ", m[x], " and the strings in ", m[x], " are:")
    count = 0
    for j in r:
      if (get_lang_abbr_from_resp(j)) == x:
         print(lang_samples[count])
      count += 1
         

def get_lang_abbr_from_resp(http_resp):
    """
    This function takes a requests object containing a response from
    detectlanguage.com, parses it, and returns the abbreviation of
    the language detected.
    """
    return http_resp.json()["data"]["detections"][0]["language"]

def get_advice():
   r = list()
   for i in range (0, 50):
      r.append(requests.get('https://api.adviceslip.com/advice'))
   advice = []
   for slip in r:
      x = slip.json()
      print(x)
      advice.append(x['slip']['advice'])
   return advice

def find_most_popular_advice(advice_list):
    """
    The `key` parameter of the `max` function takes a function to be
    applied to each element in the list. The `max` function then
    applies the key function to each element in the list with the higest
    result of the application of the key function.
    """
    return max(advice_list, key=lambda x: advice_list.count(x))



if __name__ == "__main__":
   advice = get_advice()
   var = (find_most_popular_advice(advice))

   print("most popular advice: ", var)
   seen = {}
   for i in advice:
      if i not in seen:
         seen[i] = 1
      else:
         seen[i] += 1
   
   print("times it appears: ", seen[var])

   post_language_samples()
