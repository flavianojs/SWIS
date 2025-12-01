import sys

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
         
         newline = ' %3i'*5%tuple( map(int, data_str[0:5]) ) + '%16.8f'*3%tuple( map(float_rescale, data_str[5:8]) ) + '%16.8f'%( value * float(data_str[8]) ) + '\n'
         target.write( newline )

   target.close()

print("Scale pair interactions:")

if len(sys.argv) == 1:
   #Default inputs

   # scale_input='inputfiles/pair_Mn5Ge3_rad4.00.txt'
   scale_input='inputfiles/pair_Mn5Ge3_rad0.63.txt'
   scale_input='inputfiles/pair_Mn5Ge3_rad1.20.txt'
   scale_output='inputfiles/pair_Mn5Ge3_temp.txt'

elif len(sys.argv)==3 :
   scale_input=sys.argv[1]
   scale_output=sys.argv[2]

else:
   print("Too few or too many inputs were passed to the script.")
   print("Run it with: python scale_pair.py pair_XXX.txt pair_XXX_temp.txt")
   print("or with defatult inputs: python scale_pair.py")
   exit()

print(' scale_input  =', scale_input)
print(' scale_output =', scale_output)

rescaling_interactions(scale_input, scale_output, scalling_factor )
