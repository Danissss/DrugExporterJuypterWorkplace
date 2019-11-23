import os, sys
import csv
import random





def SplitDataSet(csvFile):

	print("=========================================================================")
	print(csvFile)

	training_file = csvFile.replace(".csv","_training.csv")
	testing_file  = csvFile.replace(".csv","_testing.csv")

	training_list = []
	testing_list  = []

	training_writer = csv.writer(open(training_file, 'w', newline=''),quoting=csv.QUOTE_ALL)
	testing_writer = csv.writer(open(testing_file, 'w', newline=''),quoting=csv.QUOTE_ALL)


	attribute_name  = GetAttributeNameFake(1227)
	training_writer.writerow(attribute_name)
	testing_writer.writerow(attribute_name)

	descriptor_name = []

	# descriptorValue = []
	# save everything in list 
	# count positive and negative instance number
	total_instance = []
	negative = 0
	positive = 0
	length_descriptor = 0
	with open(csvFile) as fd:
		rd = csv.reader(fd, delimiter=",")
		descriptor_name += next(rd)
		length_descriptor = len(descriptor_name)
		for row in rd:
			total_instance.append(row)
			if "Non" in row[length_descriptor-1]:
				negative += 1
			else:
				positive += 1

	select_positive = round(positive * 0.2)
	select_negative = round(negative * 0.2)
	print("select_positive => {0}; select_negative => {1};".format(select_positive,select_negative))


	while select_negative != 0:
		random_int = random.randint(0, len(total_instance)-1)
		instance = total_instance[random_int]
		if "Non" in instance[length_descriptor-1]:
			testing_list.append(instance)
			total_instance.remove(instance)
			select_negative -= 1



	while select_positive != 0:
		random_int = random.randint(0, len(total_instance)-1)
		instance = total_instance[random_int]
		if "Non" not in instance[length_descriptor-1]:
			testing_list.append(instance)
			total_instance.remove(instance)
			select_positive -= 1

	for i in total_instance:
		training_writer.writerow(i)

	for i in testing_list:
		testing_writer.writerow(i)



def GetAttributeNameFake(size):
	attribute_name = []
	for i in range(size):
		attribute_name.append("Attribute_{0}".format(i))
	attribute_name.append("Class")
	return attribute_name






def main():
	SplitDataSet(r'C:\Users\Danis\Desktop\GetDescriptorsBindingData\result\MDR1_binding_dataset_3DFile_3D_descriptor_value.csv')
	SplitDataSet(r'C:\Users\Danis\Desktop\GetDescriptorsBindingData\result\BCRP_binding_dataset_3DFile_3D_descriptor_value.csv')
	SplitDataSet(r'C:\Users\Danis\Desktop\GetDescriptorsBindingData\result\MRP1_binding_dataset_3DFile_3D_descriptor_value.csv')
	SplitDataSet(r'C:\Users\Danis\Desktop\GetDescriptorsBindingData\result\MRP2_binding_dataset_3DFile_3D_descriptor_value.csv')















if __name__ == '__main__':
	main()