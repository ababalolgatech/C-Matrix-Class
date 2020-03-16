#ifndef _CONFIGVARIO_H_
#define _CONFIGVARIO_H_

#include<vector>
#include<map>
#include<iostream>
#include <algorithm>
#include <string>
#include <fstream>
#include <sstream>

class convert
{
public:
	template <typename T>
	static std::string T_to_string(T const &val)
	{
		std::ostringstream ostr;
		ostr << val;

		return ostr.str();
	}

	template <typename T>
	static T string_to_T(std::string const &val)
	{
		std::istringstream istr(val);
		T returnVal;
		if (!(istr >> returnVal))
			exitWithError("CFG: Not a valid " + (std::string)typeid(T).name() + " received!\n");

		return returnVal;
	}

	template <>
	static std::string string_to_T(std::string const &val)
	{
		return val;
	}

	template<typename T>
	static T StringToNumber(const std::string& numberAsString)
	{
		T valor;

		std::stringstream stream(numberAsString);
		stream >> valor;
		if (stream.fail()) {
			std::runtime_error e(numberAsString);
			throw e;
		}
		return valor;
	}
};

/*
void exitWithError(const std::string &error)
{
	std::cout << error;
	std::cin.ignore();
	std::cin.get();

	exit(EXIT_FAILURE);
}
*/



class MyMap : public std::map<std::string, std::vector<float> > {
public:
	MyMap() {}
	~MyMap() {};
	std::vector<float> numvector;





	std::vector<float> find_ret(std::string str) {
		std::map<std::string, std::vector<float> >::iterator it;

		it = std::map<std::string, std::vector<float> >::find(str);
		std::vector<float> intvec = it->second;
		std::vector<float> rec;
		for (unsigned int i = 0; i < intvec.size(); i++)
		{
			rec.push_back(intvec[i]);
		}

		return rec;
	}

	/*
	void insert(std::string str, float number){
		numvector.push_back(number);
		std::map<std::string, std::vector<float> >::insert(
		std::make_pair(str, numvector));
	}
	*/
};


class Config_vario : public  MyMap {
public:

	std::string  fp, line, key;
	std::string value;
	std::vector<int> vecpos;
	float actual_value;
	MyMap contents;
	size_t  pos;
	std::string keydelimiter, delimiter;
	char commentdelimiter;
	std::vector<float> recvector, tmpvec;


	int delimiterPos;

	Config_vario(std::string fp_, std::string keydelim, std::string delim, char comment_delim) {
		fp = fp_;
		keydelim = keydelimiter;
		delim = delimiter;
		commentdelimiter = comment_delim;
		read_file(fp);
	};
	Config_vario(std::string fp_) {
		fp = fp_;
		delimiter = ";";
		keydelimiter = "=";
		commentdelimiter = '%';

		read_file(fp);
	};

	void read_file(std::string fp) {
		std::ifstream cFile(fp);
		if (cFile.is_open())
		{
			while (getline(cFile, line))
			{
				line.erase(std::remove_if(line.begin(), line.end(), isspace), line.end());
				if (line[0] == commentdelimiter || line.empty())
					continue;

				if (line.find(delimiter) == std::string::npos) {
					getkey_single();
					inserter(key);
				}
				else
				{
					delimiterPos = line.find(keydelimiter);
					key = line.substr(0, delimiterPos);
					line.erase(0, delimiterPos + delimiter.length());

					while ((pos = line.find(delimiter)) != std::string::npos)
					{
						value = line.substr(0, pos);
			auto actual_value = convert::StringToNumber<float>(value);
						//contents.insert(key, actual_value);
						tmpvec.push_back(actual_value);
						line.erase(0, pos + delimiter.length());
					}
					inserter(key);

				}
			}
		}
		else { std::cerr << "Couldn't open config file for reading.\n"; }
	}

	void getkey_single() {
		delimiterPos = line.find(keydelimiter)     ;
		         key = line.substr(0, delimiterPos);
		       value = line.substr(delimiterPos + 1);
		tmpvec.push_back(convert::StringToNumber<float>(value));
	}
	std::vector<float> return_value(std::string key) {
		recvector.clear(); // to clear previous 
		std::map<std::string, std::vector<float> >::iterator it;
		it = contents.find(key);
		std::vector<float> intvec = it->second;
		for (unsigned int i = 0; i < intvec.size(); i++)
		{
			recvector.push_back(intvec[i]);
		}
		return recvector;

	}

	void inserter(std::string) {
		//tmpvec.push_back(number);
		contents.insert(std::make_pair(key, tmpvec));
		tmpvec.clear(); // clears it for the next key

	}

	// https://ubuntuforums.org/archive/index.php/t-822175.html
};

// if key , not found throw a  better error not assertion error
#endif