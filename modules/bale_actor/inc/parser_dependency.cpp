#include <iostream>
#include <regex>
#include <vector>
#include <fstream>
#include <string>
#include <algorithm>
#include <tuple>
#include <list>
using namespace std;

int main(int argc, char* argv[])
{
    if (argc != 2) {printf("Usage: ./parser <file_name> ($TARGET only, without .cpp)\n"); return 1;}
    string TARGET = (string)argv[1];
    string file_name = TARGET + ".cpp";
    string temp_file_name = file_name;

    auto isEmptyOrBlank = [](const std::string &s) {return s.find_first_not_of(" \t") == std::string::npos;};

    // create parsed file that is exact replica of original
    ifstream in_file(file_name);
    int pos_out = file_name.find(".cpp");
    string out_name = file_name.insert(pos_out, "_parsed_dep");
    ofstream out_file(out_name);
    string out_line;
    if (in_file && out_file) {
        while(getline(in_file, out_line)){
            out_file << out_line << "\n";
        }
    }
    out_file.close();
    in_file.close();

    /**********************************************/
    /* PHASE 1 */
    /**********************************************/
    
    vector<string> class_names;
    vector<int> class_line_nums;
    vector<int> selector_line_nums_end;

    string selector_line;
    fstream selector_file;
    int selector_line_num = 0;
    selector_file.open(out_name,ios::in);
    if (selector_file.is_open()){
        while(getline(selector_file, selector_line)){
            selector_line_num++;
            auto const regex_class = regex("class[a-zA-Z0-9 ]+:");
            smatch m_class; 
            regex_search(selector_line, m_class, regex_class); 
            for (string x: m_class) {
                auto const regex_class_name_pre = regex("[a-zA-Z0-9]+[ ]*:");
                smatch m_class_name_pre;
                regex_search(x, m_class_name_pre, regex_class_name_pre);
                for (string xx: m_class_name_pre) {
                    auto const regex_class_name_final = regex("[a-zA-Z0-9]+");
                    smatch m_class_name_final;
                    regex_search(xx, m_class_name_final, regex_class_name_final);
                    for (string name_final: m_class_name_final) {
                        class_names.push_back(name_final);
                        if(!name_final.empty()) class_line_nums.push_back(selector_line_num);
                    }
                }
            }
        }
        class_names.erase(remove_if(class_names.begin(), class_names.end(), isEmptyOrBlank), class_names.end());

        selector_file.close();
    }
    class_line_nums.push_back(selector_line_num);

    /**********************************************/
    /* PHASE 2 */
    /**********************************************/

    int final_pos_out = temp_file_name.find(".cpp");
    string final_out_name = temp_file_name.insert(final_pos_out, "_parsed");

    int class_num = 0;
    for (string i: class_names) {
        selector_file.open(out_name,ios::in);
        if (selector_file.is_open()){
            int line_num = 0;
            while(getline(selector_file, selector_line)){
                line_num++;
                if ((line_num >= class_line_nums[class_num]) && (line_num < class_line_nums[class_num+1])) {
                    if (selector_line.find(i + "(") != string::npos) {
                        fstream nested_file;
                        string nested_line;
                        int nested_line_num = 0;
                        int balanced_parenthesis = 0;
                        bool brackets_first = false;
                        nested_file.open(out_name,ios::in);
                        if (nested_file.is_open()){
                            while(getline(nested_file, nested_line)) {
                                nested_line_num++;
                                if ((nested_line_num >= line_num) & (nested_line_num < class_line_nums[class_num+1])) {
                                    auto const regex_brackets_left =  regex("[{]+");
                                    auto const regex_brackets_right =  regex("[}]+");
                                    smatch m_brackets_left;
                                    smatch m_brackets_right;
                                    regex_search(nested_line, m_brackets_left, regex_brackets_left);
                                    regex_search(nested_line, m_brackets_right, regex_brackets_right);

                                    for (string t: m_brackets_left) balanced_parenthesis += t.size();
                                    
                                    if (balanced_parenthesis > 0) brackets_first = true;

                                    for (string t: m_brackets_right) balanced_parenthesis -= t.size();
                                    
                                    if (balanced_parenthesis==0 & brackets_first) {
                                        selector_line_nums_end.push_back(nested_line_num);
                                        break;
                                    }
                                }
                            }
                        }
                        nested_file.close();
                    }
                }
            }
            selector_file.close();
        }

        vector<string> mb_ids;
        vector<string> mb_process_funcs;
        string phase2_line;
        fstream phase2_file;
        phase2_file.open(out_name,ios::in);
        int phase2_line_num = 0;
        if (phase2_file.is_open()){
            while(getline(phase2_file, phase2_line)){
                phase2_line_num++;
                if ((phase2_line_num >= class_line_nums[class_num]) && (phase2_line_num < class_line_nums[class_num+1])) {
                    // get mb id and proc functions
                    auto const regex_mb = regex("mb\\[[a-zA-Z0-9]+\\].process");
                    auto const regex_process = regex("this->[a-zA-Z0-9_]+\\(");
                    smatch m_mb; 
                    smatch m_process; 
                    regex_search(phase2_line, m_mb, regex_mb); 
                    regex_search(phase2_line, m_process, regex_process); 

                    for (string x : m_mb) {
                        auto const regex_mb_pre_final =  regex("[a-zA-Z0-9]*\\]");
                        smatch m_mb_pre_final;
                        regex_search(x, m_mb_pre_final, regex_mb_pre_final);
                        for (string xx: m_mb_pre_final) {
                            auto const regex_mb_final =  regex("[a-zA-Z0-9]+");
                            smatch m_mb_final;
                            regex_search(xx, m_mb_final, regex_mb_final);
                            for (string mb_final: m_mb_final)
                                if (!mb_final.empty()) mb_ids.push_back(mb_final);
                        }
                    }
                    mb_ids.erase(remove_if(mb_ids.begin(), mb_ids.end(), isEmptyOrBlank), mb_ids.end());

                    for (string x : m_process) {
                        auto const regex_process_pre_final =  regex("[a-zA-Z0-9_]*\\(");
                        smatch m_process_pre_final;
                        regex_search(x, m_process_pre_final, regex_process_pre_final);
                        for ( string xx: m_process_pre_final) {
                            auto const regex_process_final =  regex("[a-zA-Z0-9_]+");
                            smatch m_process_final;
                            regex_search(xx, m_process_final, regex_process_final);
                            for ( string process_final: m_process_final)
                                if (!process_final.empty()) mb_process_funcs.push_back(process_final);
                        }
                    }
                    mb_process_funcs.erase(remove_if(mb_process_funcs.begin(), mb_process_funcs.end(), isEmptyOrBlank), mb_process_funcs.end());
                }
            }
            phase2_file.close();
        }

        vector<string> mb_dependencies(mb_ids.size());
        // dependencies case 1: can be within mb[].process {} definition
        phase2_file.open(out_name,ios::in);
        phase2_line_num = 0;
        int index_m = 0;
        if (phase2_file.is_open()){
            while(getline(phase2_file, phase2_line)){
                phase2_line_num++;
                if ((phase2_line_num >= class_line_nums[class_num]) && (phase2_line_num < class_line_nums[class_num+1])) {
                    auto const regex_mb = regex("mb\\[[a-zA-Z0-9]+\\].process");
                    smatch m_mb; 
                    regex_search(phase2_line, m_mb, regex_mb); 

                    for (string x : m_mb) {
                        auto const regex_mb_pre_final =  regex("[a-zA-Z0-9]*\\]");
                        smatch m_mb_pre_final;
                        regex_search(x, m_mb_pre_final, regex_mb_pre_final);
                        for (string xx: m_mb_pre_final) {
                            auto const regex_mb_final =  regex("[a-zA-Z0-9]+");
                            smatch m_mb_final;
                            regex_search(xx, m_mb_final, regex_mb_final);
                            for (string mb_final: m_mb_final) {
                                if (!mb_final.empty()) {
                                    int first_entry = true;
                                    fstream nested_file;
                                    string nested_line;
                                    int nested_line_num = 0;
                                    int balanced_parenthesis = 0;
                                    bool brackets_first = true;
                                    nested_file.open(out_name,ios::in);
                                    if (nested_file.is_open()){
                                        while(getline(nested_file, nested_line)) {
                                            nested_line_num++;
                                            if ((nested_line_num >= phase2_line_num) && (nested_line_num < class_line_nums[class_num+1])) {
                                                if (balanced_parenthesis==0 & brackets_first) brackets_first = false;
                                                auto const regex_brackets_left =  regex("[{]+");
                                                auto const regex_brackets_right =  regex("[}]+");
                                                smatch m_brackets_left;
                                                smatch m_brackets_right;
                                                regex_search(nested_line, m_brackets_left, regex_brackets_left);
                                                regex_search(nested_line, m_brackets_right, regex_brackets_right);

                                                for (string t: m_brackets_left) balanced_parenthesis += t.size();

                                                for (string t: m_brackets_right) balanced_parenthesis -= t.size();
                                                
                                                if (balanced_parenthesis==0 & !brackets_first) break;

                                                //search for send to other mailbox
                                                auto const regex_send =  regex("send\\([a-zA-Z0-9,]+");
                                                smatch m_send;
                                                regex_search(nested_line, m_send, regex_send);
                                                for (string send_final: m_send) {
                                                    if(!send_final.empty() && first_entry) {mb_dependencies[index_m] = mb_ids[index_m]+send_final; first_entry = false;}
                                                    else if(!send_final.empty() && !first_entry) {mb_dependencies[index_m] = mb_dependencies[index_m] + "/" + mb_ids[index_m]+send_final;}
                                                }
                                            }
                                        }
                                    }
                                    nested_file.close();
                                    index_m++;
                                }

                            }
                        }
                    }

                }
            }
            phase2_file.close();
        }

        // dependencies case 2: can be within void process function
        selector_file.open(out_name,ios::in);
        if (selector_file.is_open()){
            // get dependency in proc funcs
            int line_num = 0;
            while(getline(selector_file, selector_line)){
                line_num++;
                int index = 0;
                if ((line_num >= class_line_nums[class_num]) && (line_num < class_line_nums[class_num+1])) {
                for (string proc_func : mb_process_funcs) {
                    int first_entry = true;
                    if (selector_line.find("void " + proc_func) != string::npos) {
                        fstream nested_file;
                        string nested_line;
                        int nested_line_num = 0;
                        int balanced_parenthesis = 0;
                        bool brackets_first = true;
                        nested_file.open(out_name,ios::in);
                        if (nested_file.is_open()){
                            while(getline(nested_file, nested_line)) {
                                nested_line_num++;
                                if ((nested_line_num >= line_num) && (nested_line_num < class_line_nums[class_num+1])) {
                                    if (balanced_parenthesis==0 & brackets_first) brackets_first = false;
                                    auto const regex_brackets_left =  regex("[{]+");
                                    auto const regex_brackets_right =  regex("[}]+");
                                    smatch m_brackets_left;
                                    smatch m_brackets_right;
                                    regex_search(nested_line, m_brackets_left, regex_brackets_left);
                                    regex_search(nested_line, m_brackets_right, regex_brackets_right);

                                    for (string t: m_brackets_left) balanced_parenthesis += t.size();

                                    for (string t: m_brackets_right) balanced_parenthesis -= t.size();
                                    
                                    if (balanced_parenthesis==0 & !brackets_first) break;

                                    //search for send to other mailbox
                                    auto const regex_send =  regex("send\\([a-zA-Z0-9,]+");
                                    smatch m_send;
                                    regex_search(nested_line, m_send, regex_send);
                                    for (string send_final: m_send) {
                                        if(!send_final.empty() && first_entry) {mb_dependencies[index] = mb_ids[index]+send_final; first_entry = false;}
                                        else if(!send_final.empty() && !first_entry) {mb_dependencies[index] = mb_dependencies[index] + "/" + mb_ids[index]+send_final;}
                                    }
                                }
                            }
                        }
                        nested_file.close();
                    }
                    index++;
                }
                }
            }
            selector_file.close();
        }

        vector<list<string>> mb_successors;
        vector<list<string>> mb_predecessors(mb_ids.size());
        int index_ids = 0;
        for (string mb_id_f: mb_ids){
            list<string> temp_list = {};

            vector<string> mb_dependencies_split;
            string dep_name_line = mb_dependencies[index_ids];
            string delimiter = "/";
            int pos = 0;
            string token;
            while ((pos = dep_name_line.find(delimiter)) != string::npos) {
                token = dep_name_line.substr(0, pos);
                mb_dependencies_split.push_back(token);
                dep_name_line.erase(0, pos + delimiter.length());
            }
            mb_dependencies_split.push_back(dep_name_line);


            for (string mb_id_f_send: mb_dependencies_split) {
                string delimiter = "send(";
                string first_s = mb_id_f_send.substr(0, mb_id_f_send.find(delimiter));   
                mb_id_f_send.erase(0, mb_id_f_send.find(delimiter) + delimiter.length());  

                if(first_s == mb_id_f) {
                    delimiter = ",";
                    string second_s = mb_id_f_send.substr(0, mb_id_f_send.find(delimiter)); 
                    temp_list.push_front(second_s);
                    auto it = find(mb_ids.begin(), mb_ids.end(), second_s);
                    int index_dep = it - mb_ids.begin();
                    mb_predecessors[index_dep].push_front(first_s);
                }
            }
            mb_successors.push_back(temp_list);
            index_ids++;
        }  

        for (int i = 0; i < mb_successors.size(); i++) {
            mb_successors[i].sort(); mb_successors[i].unique();
            mb_predecessors[i].sort(); mb_predecessors[i].unique();
        }

        ifstream final_in_file(out_name);
        ofstream final_out_file(final_out_name);
        string final_out_line;
        int final_line_num = 0;
        if (final_in_file && final_out_file) {
            while(getline(final_in_file, final_out_line)){
                final_line_num++;
                final_out_file << final_out_line << "\n";
                if (final_line_num == (selector_line_nums_end[class_num]-1)) {
                    int mb_final_id_index = 0;
                    for (string mb_final_id: mb_ids) {
                        string final_predecessors;
                        string final_successors;
                        bool first_predecessor = true;
                        bool first_successor = true;

                        for (string u: mb_predecessors[mb_final_id_index]) {
                            if (first_predecessor) {final_predecessors = u; first_predecessor = false;}
                            else {final_predecessors = final_predecessors + ", " + u;}
                        }

                        for (string u: mb_successors[mb_final_id_index]) {
                            if (first_successor) {final_successors = u; first_successor = false;}
                            else {final_successors = final_successors + ", " + u;}
                        }

                        string final_line = "mb[" + mb_final_id + "].add_dep_mailboxes({" + final_predecessors + "}, {" + final_successors + "});";
                        final_out_file << final_line << "\n";

                        mb_final_id_index++;
                    }
                }
            }
        }
        final_out_file.close();
        final_in_file.close();

        //copy all changes back
        ifstream in_file(final_out_name);
        ofstream out_file(out_name);
        string out_line;
        if (in_file && out_file) {
            while(getline(in_file, out_line)){
                out_file << out_line << "\n";
            }
        }
        out_file.close();
        in_file.close();

        class_num++;
    }

    return 0; 
}