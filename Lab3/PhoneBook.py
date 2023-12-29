####################################################################################################################
###                                    Algorithms for Bioinformatics
###                                 ***  class Phone Book              ***
###
### Test number: 4      Class Number: 3         Date:   2022
###
### Group
### Student: ....               Number:...
### Student: ....               Number:...
###
####################################################################################################################
### Complete the code below for the object PhoneBook
### In main give example on how to create, update, insert and use object PhoneBook
### Explain in comments how the data will be organized

class PhoneBook:
    ''' Implements a Phone Book '''


    def __init__(self,myName,myNumber):
        ''' initializes phone book with appropriate data structure '''
        self.dic={}
        self.myName=myName
        self.myNumber=myNumber

    def add_phone(self, name, number):
        if name in self.dic.keys():
            if type(self.dic[name])=="list":
                self.dic[name]=self.dic[name]+[number]
            else:
                self.dic[name]=[self.dic[name],number]
        else:
            self.dic[name]=number

    def search_by_name(self,name):
        if name in self.dic.keys():
            return self.dic[name]
        else:
            return "No number with that name associated!"

    def search_by_phoneNumber(self,number):
        if number in self.dic.items():
            return self.dic.keys()[self.dic.values().index(number)]
        else:
            return "No name with that number associated!"

    def print_phoneBook(self):
        print("Phone Book    Owner:",self.myName,"|OwnerNumber:",self.myNumber)
        print("{:<15}  {:<15}".format("Name","Number"))
        for i in self.dic.keys():
            print("{:<15}: {:<15}".format(i,str(self.dic[i])))

    def copy_phoneBook(self):
        return self.dic


if __name__ == "__main__":
    #test code here
    myName="Y"
    myNumber="123"
    pB=PhoneBook(myName,myNumber)
    pB.add_phone("abc","123124")
    pB.add_phone("abc","123")
    pB.add_phone("asdfsdf","12435245")
    print(pB.search_by_name("abc"))
    print(pB.search_by_name("asdfasdfasdf"))
    print(pB.search_by_phoneNumber("123"))
    print(pB.search_by_phoneNumber("3333"))
    print(pB.copy_phoneBook())
    pB.print_phoneBook()
