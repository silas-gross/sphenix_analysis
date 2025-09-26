#include <fstream>
#include <string>
#include <sstream>

//just grab the number of events and scale the cross section

int GetFilteredCrossSection(std::string log_file, std::string filterlog, std::string filter_cross_section_file)
{
	fstream inital_file (log_file);
	fstream filtered_file (filterlog);
	fstream output_file (fiter_cross_section_file);
	std::string inital_cross_section, surviving_fraction;
	std::string log_line, filter_line;
	while(std::get_line(inital_file, log_line))
	{
		if(log_line.find("generated events): ") != std::string::npos)
		{
			std::stringstream cross_sect (log_line);
			std::string cs_s;
			while(std::get_line(cross_sect, cs_s, ' '))
			{
				if(cs_s.find("e+") !=std::string::npos || cs_s.find("e-") != std::string::npos)
				{
					initial_cross_section=cs_s;
				}
			}
		}
	}
	initial_file.close();
	while(std::get_line(filtered_file, filter_line))
	{
		

