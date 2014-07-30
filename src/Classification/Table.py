'''
Created on Sep 16, 2010

@author: jhooge
'''
"""
Parse a file and put it into a table structure
"""

def evaluate(word):
    try:
        return eval(word)
    except:
        return word

class Table(list):
    '''
    classdocs
    '''

    def load(self, lines, format = None, sep=None):
        """
        converts a list of text lines into a table
        """
        for line in lines:
            if format is None:
                words = line.split(sep)
            else:
                words = [line[start:end].strip() for (start, end) in format]
            words = map(evaluate, words)
            self.append(words)

    def read(self, filename, format = None, sep=None):
        """
        reads file and puts the content into a table
        """
        import os
        file = open(os.path.expanduser(filename))
        lines = file.readlines()
        file.close()
        self.load(lines, format, sep)

    def rows(self, indices):
        return Table([self[index] for index in indices])

    def columns(self, indices):
        return Table([[self[i][j] for j in indices]
                      for i in range(len(self))])

    def del_rows(self, indices):
        indices.sort()
        indices.reverse()
        for i in indices: del self[i]

    def del_columns(self, indices):
        for i in range(len(self)):
            for j in indices:
                del self[i][j]

    def write(self, filename, separator = ' '):
        import os

        file = open(os.path.expanduser(filename),'w')
        [file.write(separator.join(map(str, self[i])) + '\n')
         for i in range(len(self))]
        file.close()
