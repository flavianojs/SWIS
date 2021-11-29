scalling_factor = 1.0

def float_rescale(x) :
   return float(x) * scalling_factor

def rescaling_interactions(filename, filename_temp, value):
   # openning and reading the file
   file = open(filename, 'r')
   lines=file.readlines()
   file.close()

   # reopening the file in the writing mode
   target = open(filename_temp, 'w')

   # looping through the files looking for keywords
   for line_index, line in enumerate(lines) :
      
      if line_index == 0 :
         target.write( line )
      
      else :
         data_str=line.split()
         
         newline = '%3i'*5%tuple( map(int, data_str[0:5]) ) + '%16.8f'*3%tuple( map(float_rescale, data_str[5:8]) ) + '%16.8f'%( value * float(data_str[8]) ) + '\n'
         target.write( newline )

   target.close()

rescaling_interactions('inputfiles/pair_MnSi_rad4.5noncol.txt', 'inputfiles/pair_MnSi_rad4.5noncol_temp.txt', scalling_factor )