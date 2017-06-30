newlines=[]
with open('em_close.itp','r') as f:
	for line in f:
		try:
			num1 =int(line.partition(' ')[0])
		except ValueError:
			continue
		num2 = num1-9173
		newlines.append(line.replace(str(num1),str(num2)))

with open('em_close_renumbered.itp','w') as f:
	f.write('[ position_restraints ]\n')
	for line in newlines:
		f.write(line)
