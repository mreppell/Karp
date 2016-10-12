#include "common_k.h"

using namespace std;

std::string pretty_num(int num) {
  return pretty_num(static_cast<size_t>(num));
}

std::string pretty_num(unsigned int num) {
  return pretty_num(static_cast<size_t>(num));
}

std::string pretty_num(size_t num) {
  auto s = std::to_string(num);
  auto ret = std::string("");

  if (s.size() <= 3) {
    return s;
  }

  int remainder = s.size() % 3;
  if (remainder == 0) {
    remainder = 3;
  }

  size_t start_pos = 0;
  while (start_pos + remainder < s.size() - 1) {
    ret += s.substr(start_pos, remainder) + ",";
    start_pos += remainder;
    remainder = 3;
  }

  ret += s.substr(start_pos, 3);

  return ret;
}

void split(const string& str, const string& delimiters, vector<string>& tokens) {

  if (!str.empty()) {
    size_t d_w = delimiters.length();
    size_t s_w = str.length();

    //string::size_type lastPos = str.find(delimiters);
    string::size_type lastPos = 0;  
    string::size_type pos = str.find(delimiters);
    if (pos==string::npos) {
      tokens.push_back(str);
    } else {
      while (string::npos != pos) {
	tokens.push_back(str.substr(lastPos,pos-lastPos));
	lastPos=pos+d_w;
	pos = str.find(delimiters,lastPos);
      }
      if(lastPos < s_w) {
	tokens.push_back(str.substr(lastPos,s_w));
      }
    }
  }
}
