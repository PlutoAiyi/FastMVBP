#include <algorithm>
#include <cstdio>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>
#include <omp.h>
using std::max;
using std::min;
using std::unordered_map;
#include "Framework.hpp"
#include "Graph.hpp"
extern udi Framework::k;
extern udi Framework::top_k;
void transform(std::string input_filename, std::string output_filename,int all_size,int left_part_size,int edges) {
    std::ifstream input_file(input_filename);
    std::ofstream output_file(output_filename);

    if (!input_file.is_open()) {
        std::cerr << "Error: Unable to open input file: " << input_filename << std::endl;
        return;
    }

    if (!output_file.is_open()) {
        std::cerr << "Error: Unable to open output file: " << output_filename << std::endl;
        return;
    }
    std::string line;
    std::vector<std::vector<int>> output;
    output.resize(all_size);
    for(int i=0;i<all_size;i++){
        output[i].push_back(i);
    }
    int last_node=-1,last_left=-1;
    while (std::getline(input_file, line)) {
        if (line.empty() || line[0] == '%') {
            continue;
        }

        std::istringstream iss(line);
        int left_node, right_node;
        if (!(iss >> left_node >> right_node)) {
            std::cerr << "Error: Invalid input format." << std::endl;
            return;
        }
        
        // Adjust node indices to start from 0
        left_node -= 1;
        right_node -=1 ;
        right_node +=left_part_size;
        if(last_node==right_node && last_left==left_node) continue;
        output[left_node].push_back(right_node);
        output[right_node].push_back(left_node);
        last_node=right_node;
        last_left=left_node;
    }
    for (int i = 0; i < all_size; ++i) {
        std::sort(output[i].begin()+1, output[i].end());
        output[i].erase(std::unique(output[i].begin()+1, output[i].end()), output[i].end());
    }
    std::ofstream outputFile(output_filename);

    if (outputFile.is_open()) {
        outputFile<< all_size <<" "<<left_part_size<<" "<<edges<< std::endl;
        for (int i = 0; i < all_size; ++i) {
            for (size_t j = 0; j < output[i].size()-1; ++j) {
                outputFile << output[i][j] << " ";
            }

            outputFile<< output[i][output[i].size()-1] << std::endl;
        }

        // 关闭文件
        outputFile.close();
        std::cout << "Data has been written to output.txt." << std::endl;
    } else 
        std::cerr << "Unable to open the file." << std::endl;
    }
int main(int argc, char* argv[]) {
    int k = -1;
    int top_k = -1;
    std::string graphPath = "";
    bool gFlag = false; 
    int theta_L = -1;
    int theta_R = -1;
    int all_size=-1;
    int left_size=-1;
    int edges=-1;
    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        if (arg == "-k" && i + 1 < argc) {
            k = std::stoi(argv[++i]);
        } else if (arg == "-K" && i + 1 < argc) {
            top_k = std::stoi(argv[++i]);
        } else if (arg == "-f" && i + 1 < argc) {
            graphPath = argv[++i];
        } else if (arg == "-l" && i + 1 < argc) {
            theta_L = std::stoi(argv[++i]);
        } else if (arg == "-r" && i + 1 < argc) {
            theta_R = std::stoi(argv[++i]);
        } else if (arg == "-g" && i + 3 < argc) {
            gFlag = true;
            all_size = std::stoi(argv[++i]);
            left_size = std::stoi(argv[++i]);
            edges = std::stoi(argv[++i]);
        } else {
            std::cerr << "Unknown argument or missing value: " << arg << std::endl;
            return 1;
        }
    }
    if (gFlag) {
        std::cout << "Entered the graph processing branch" << std::endl;
        if (graphPath.empty()) {
            std::cerr << "Error: -f option is required when using -g option" << std::endl;
            return 1;
        }
        std::string new_filepath = graphPath;
        size_t last_dot_position = new_filepath.find_last_of('.');
        if (last_dot_position != std::string::npos) {
            new_filepath = new_filepath.substr(0, last_dot_position);
        }
        new_filepath += ".g";
        std::cout << "Graph size: " << all_size << ", left size: " << left_size << ", edges:" << edges << std::endl;
        transform(graphPath,new_filepath,all_size,left_size,edges);
        std::cout<<"transform success"<<std::endl;
    }else{
        if (k == -1 || top_k == -1 || graphPath.empty()) {
            std::cerr << "Missing required arguments" << std::endl;
            return 1;
        }
        if (theta_L != -1 && theta_R != -1 && (theta_L < 2 * k + 1 || theta_R < 2 * k + 1)){
            std::cerr << "constraints have to be larger than 2k+1 " << std::endl;
            return 1;
        }
        if (theta_L == -1) theta_L = 2 * k + 1;
        if (theta_R == -1) theta_R = 2 * k + 1;
        Framework::k = k;
        Framework::top_k = top_k;
        std::cout<<"start reading "<<graphPath<<std::endl;
        Graph g(graphPath);
        std::cout<<"reading finish"<<std::endl;
        std::cout<<"Left vertexs:"<<g.L_size<<",Right vertexs:"<<g.vertex_tot-g.L_size<<",Edges:"<<g.edge_tot<<std::endl;
        time_t start = clock();
            Framework::IB(g, theta_L, theta_R); 
        std::cout << "tot time = " << (double)(clock() - start) / CLOCKS_PER_SEC << " seconds"<<std::endl;

        udi cnt = 1;
        std::cout << "tot cnt = " << Framework::pq_k.size() << std::endl;
        while (!Framework::pq_k.empty()){
            Framework::pq_k.top().print();
            Framework::pq_k.pop();
            cnt++;
        }
    }
    return 0;
}