#include <iostream>
#include <sys/stat.h>
#include <sys/types.h>
#include <random>
#include <chrono>
#include <cstring>
#include <set>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <memory>
#include <unordered_set>
#include <unordered_map>
#include <map>
#include <limits>
#include "noisyBinarySearch.hpp"


#define C 4 // used during RandomizedShellSort

long long num_comparisons;

std::vector<int> input, output;
int n;
int comparison_cnt=0;

std::uniform_real_distribution<double> distribution(0, 1);

unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
std::default_random_engine generator (seed);

// doubly linked list used for window-merge sort
struct Node {
  int data;
  struct Node* next;
  struct Node* prev;
};

class dll {
  public:
    int size = 0;
    struct Node* head = NULL;
    struct Node* tail = NULL;
    void insertEnd(int);
    void deleteNode(struct Node*);
};

void dll::insertEnd(int data) {
  struct Node* newNode = new Node;
  size++;

  newNode->data = data;
  newNode->next = NULL;

  if (head == NULL) {
    newNode->prev = NULL;
    head = newNode;
    tail = newNode;
    return;
  }

  tail->next = newNode;
  newNode->prev = tail;
  tail = newNode;
}

void dll::deleteNode(struct Node* del_node) {
  if (head == NULL || del_node == NULL)
    return;

  if (head == del_node)
    head = del_node->next;

  if (del_node->next != NULL)
    del_node->next->prev = del_node->prev;

  if (del_node->prev != NULL)
    del_node->prev->next = del_node->next;

  free(del_node);
  size--;
}

int hash_table[8][256];
int number_of_hash_output_bits = 14;

void fill_hash_table(){ //number_of_output_bits should be less than 14
  for (int i = 0; i < 7 ; i++)
    for (int j = 0; j < (1 << 8); j++){
      hash_table[i][j] = rand() % (1 << number_of_hash_output_bits);
    }
}

int lookup_hash_value(long long number){
  int res = 0;
  for (int i = 0; i < 8; i++){
    res ^= hash_table[i][(number >> (8 * i)) % (1 << 8)];
  }
  return res;
}

std::vector<int> random_sequence_generator(std::vector<int> &v){
  for(int i=v.size()-1;i>0;i--){
    int pos = rand() % (i+1);
    std::swap(v[i], v[pos]);
  }
  return v;
}

bool pass_through_noise(int x, int y, double p){
  if(x>y)
    std::swap(x, y);
  long long conctated_number = (x * (1ll << 32)) | y;
  double hash_value = lookup_hash_value(conctated_number);
  double scaled_p = p * (1 << number_of_hash_output_bits);
  return hash_value >= scaled_p;
}

bool handle_inversion_with_noise(std::vector<int>& v, int i, int j, double p){
  bool is_truthful = pass_through_noise(v[i], v[j], p);
  num_comparisons++;
  comparison_cnt++;
  if((is_truthful && v[i]>v[j]) || (!is_truthful && v[i]<v[j])){
    std::swap(v[i], v[j]);
    return true;
  }
  return false;
}
bool is_inversion_with_noise_val(int x, int y, double p){
  num_comparisons++;
  bool is_truthful = pass_through_noise(x, y, p);
  if((is_truthful && x>y) || (!is_truthful && x<y)){
    return true;
  }
  return false;
}

bool is_inversion_with_noise(std::vector<int>& v, int i, int j, double p){
  num_comparisons++;
  bool is_truthful = pass_through_noise(v[i], v[j], p);
  comparison_cnt++;
  if((is_truthful && v[i]>v[j]) || (!is_truthful && v[i]<v[j])){
    return true;
  }
  return false;
}

long long getMaximumDisloc(std::vector<int> &A) {
  std::vector<std::pair<int, int> > A_pairs, A_sorted(A.size());
  for (int i = 0; i<A.size(); i++) {
    A_pairs.push_back(std::make_pair(A[i], i));
  }
  partial_sort_copy(begin(A_pairs), end(A_pairs), begin(A_sorted), end(A_sorted));
  long long max_disloc = 0;
  for(long long i = 0; i < A_sorted.size(); i++){
    max_disloc = std::max(max_disloc, std::abs(A_sorted[i].second - i));
  }
  return max_disloc;
}

void insertion_sort(std::vector<int> & v, double p) {
  for(int i = 1; i < v.size(); i++){
    for(int j = i-1; j >= 0; j--){
      if(!handle_inversion_with_noise(v, j, j+1, p))
        break;
    }
  }
}

int partition (std::vector<int>& arr, int low, int high, double p) {
    int i = (low - 1);
    for (int j = low; j <= high- 1; j++) {
        if (!is_inversion_with_noise(arr, j, high, p)) {
            i++;
            std::swap(arr[i], arr[j]);
        }
    }
    std::swap(arr[i + 1], arr[high]);
    return (i + 1);
}

void quickSort(std::vector<int>& arr, int low, int high, double p) {
    if (low < high) {
        int pivot = partition(arr, low, high, p);
        quickSort(arr, low, pivot - 1, p);
        quickSort(arr, pivot + 1, high, p);
    }
}

void compareRegions(std::vector<int> &A, int s, int t, int offset, double p) {
  std::vector<int> mate(offset, 0);
  for (int count = 0; count < C; count++) {
    for (int i = 0; i < offset; i++) {
      mate[i] = i;
    }
    std::random_shuffle(std::begin(mate), std::end(mate));
    for (int i = 0; i < offset; i++) {
      if(s+i < t+mate[i])
        handle_inversion_with_noise(A, s + i, t + mate[i], p);
      if(s+i> t+mate[i])
        handle_inversion_with_noise(A, t + mate[i], s+i, p);
    }
  }
}

std::vector<int> generate_pratt(int max_size) {
  std::vector<int> products;
  int pow3 = 1;
  while (pow3 <= max_size) {
    int pow2 = pow3;
    while (pow2 <= max_size) {
      products.push_back(pow2);
      pow2 *= 2;
    }
    pow3 *= 3;
  }
  std::sort(products.begin(), products.end());
  return products;
}

int shellSortPratt(std::vector<int>& arr, double p)
{
  int sz = arr.size();
  std::vector<int> sequence = generate_pratt(sz);
  for (auto gap: sequence) {
    for (int i = gap; i < sz; i += 1) {
      int temp = arr[i];
      int j;
      for (j = i; j >= gap; j -= gap) {
        if (!is_inversion_with_noise_val(arr[j-gap], temp, p)) {
          break;
        }
        arr[j] = arr[j - gap];
      }
      arr[j] = temp;
    }
  }
  return 0;
}

void randomizedShellSortWithout2s3s(std::vector<int> &A,  double p) {
  int sz = A.size();
  for (int offset = sz/2; offset > 0; offset /= 2) {
    for (int i = 0; i < sz - offset; i += offset) { // compare-exchange up
      compareRegions(A, i, i+offset, offset, p);
    }
    for (int i = sz-offset; i >= offset; i -= offset) { // compare-exchange down
      compareRegions(A, i - offset, i, offset, p);
    }
    for (int i = 0; i < sz; i += 2*offset) { // compare odd-even regions
      compareRegions(A, i, i + offset, offset, p);
    }
    for (int i = 0; i < sz - offset; i += 2*offset) { // compare even-odd regions
      compareRegions(A, i, i + offset, offset, p);
    }
  }
}

void randomizedShellSort(std::vector<int> &A,  double p) {
  int sz = A.size();
  for (int offset = sz/2; offset > 0; offset /= 2) {
    for (int i = 0; i < sz - offset; i += offset) { // compare-exchange up
      compareRegions(A, i, i+offset, offset, p);
    }
    for (int i = sz-offset; i >= offset; i -= offset) { // compare-exchange down
      compareRegions(A, i - offset, i, offset, p);
    }
    for (int i = 0; i < sz - 3*offset; i += offset) { // compare 3 hops up
      compareRegions(A, i, i + 3*offset, offset, p);
    }
    for (int i = 0; i < sz - 2*offset; i += offset) { // compare 2 hops up
      compareRegions(A, i, i + 2*offset, offset, p);
    }
    for (int i = 0; i < sz; i += 2*offset) { // compare odd-even regions
      compareRegions(A, i, i + offset, offset, p);
    }
    for (int i = 0; i < sz - offset; i += 2*offset) { // compare even-odd region
      compareRegions(A, i, i + offset, offset, p);
    }
  }
}

// generates random number in range [l, r]
int getRandomNumber(int l, int r) {
  if (l > r) {
    std::swap(l, r);
  }
  std::random_device dev;
  std::mt19937 rng(dev());
  std::uniform_int_distribution<std::mt19937::result_type> dist6(l, r);
  return dist6(rng);
}

std::pair<std::vector<int>, std::vector<int> > getAnnealingSchedule(int n, int q, double g_scale, int c, int r) {
  std::vector<int> T, R;
  // phase 1
  for (int i = n/2; i >= q*std::pow(std::log2(n), 6); i/=2) {
    T.push_back(i);
    T.push_back(i);
    R.push_back(c);
    R.push_back(c);
  }

  // phase 2
  int glogn = std::floor(g_scale*64*std::pow(std::exp(1.0), 2)*std::log2(n))+1;
  for (int i = q*std::pow(std::log2(n), 6); i >= glogn; i/=2) {
    T.push_back(i);
    R.push_back(r);
  }
  // phase 3

  for (int i = 0; i < glogn; i++) {
    T.push_back(1);
    R.push_back(1);
  }
  return std::make_pair(T, R);
}

// T: temparature sequence, R: repetition sequence
void annealingSort(std::vector<int> &A, double p, double g_scale, double h) {
  std::pair<std::vector<int>, std::vector<int> > schedule = getAnnealingSchedule(A.size(), 1, g_scale, 10, h*std::log2(n)/std::log2(std::log2(n)));
  std::vector <int> T = schedule.first;
  std::vector <int> R = schedule.second;
  int t = T.size(), sz = A.size();
  for (int j = 0; j < t; j++) {
    for (int i = 0; i < sz-1; i++) {
      for (int k = 0; k < R[j]; k++) {
        int s = getRandomNumber(i + 1, std::min(sz - 1, i + T[j]));
        if(i<s)
          handle_inversion_with_noise(A, i, s, p);
        if(s<i)
          handle_inversion_with_noise(A, s, i, p);
      }
    }
    for (int i = sz-1; i > 0; i--) {
      for (int k = 0; k < R[j]; k++) {
        int s = getRandomNumber(std::max(0, i - T[j]), i - 1);
        if(i<s)
          handle_inversion_with_noise(A, i, s, p);
        if(s<i)
          handle_inversion_with_noise(A, s, i, p);
      }
    }
  }
}

inline bool exists (const std::string& name) {
  struct stat buffer;
  return (stat (name.c_str(), &buffer) == 0);
}


void windowSort(std::vector<int>& v, int l, int r, int d1, int d2, double p){
  std::vector<int> v_temp = std::vector<int>(v.begin() + l, v.begin() + r), v_res;
  std::vector<std::pair<int,int> > scores;
  for(int w = 2*d1; w>=2*d2;w/=2){
    scores.clear();
    for(int i=0;i<v_temp.size();i++){
      int score_w = std::max(0, i-w);
      for(int j = std::max(0, i-w); j< i; j++){
        if(!is_inversion_with_noise(v_temp, j, i, p))
          score_w++;
      }
      for(int j = std::min( i + w, int(v_temp.size())-1); j > i; j--){
        if(is_inversion_with_noise(v_temp, i, j, p))
          score_w++;
      }
      scores.push_back(std::make_pair(score_w, v_temp[i]));
    }
    sort(scores.begin(), scores.end());
    v_res.clear();
    for(int i=0;i<scores.size();i++)
      v_res.push_back(scores[i].second);
    v_temp = v_res;
  }
  for (int i = l, j = 0; i < r && j < v_res.size(); i++, j++) {
    v[i] = v_res[j];
  }
}


void windowMerge(std::vector<int> &A, int l, int m, int r, double p, double d) {
  std::vector<int> v1 = std::vector<int>(A.begin(), A.begin() + m); std::vector<int> v2 = std::vector<int>(A.begin()+m, A.begin() + r);
  int test1 = getMaximumDisloc(v1);
  int test2 = getMaximumDisloc(v2);
  dll A1, A2;
  std::vector<int> output;
  std::unordered_map<int, struct Node*> valToNode;
  for (int i = l; i<m; i++) {
    A1.insertEnd(A[i]);
  }
  for (int i = m; i<r; i++) {
    A2.insertEnd(A[i]);
  }
  int a = A.size();
  while (A1.size + A2.size > 6*d) {
    std::vector<int> S;
    std::unordered_map<int, std::pair<struct Node*, int> > valToNode;
    struct Node *curr1 = A1.head, *curr2 = A2.head;
    for (int i = 0; i < std::min<double>(3*d, A1.size); i++){
      S.push_back(curr1->data);
      valToNode[curr1->data] = std::make_pair(curr1, 1);
      curr1 = curr1->next;
    }
    for (int i = 0; i < std::min<double>(3*d, A2.size); i++){
      S.push_back(curr2->data);
      valToNode[curr2->data] = std::make_pair(curr2, 2);
      curr2 = curr2->next;
    }
    windowSort(S, 0, S.size(), 4*d, d, p);
    for (int i = 0; i < d; i++) {
      output.push_back(S[i]);
      auto pair = valToNode[S[i]];
      if (pair.second == 1){
        A1.deleteNode(pair.first);
      } else {
        A2.deleteNode(pair.first);
      }
    }
  }
  std::vector<int> S;
  struct Node *curr1 = A1.head, *curr2 = A2.head;
  for (int i = 0; i < A1.size; i++){
    S.push_back(curr1->data);
    curr1 = curr1->next;
  }
  for (int i = 0; i < A2.size; i++){
    S.push_back(curr2->data);
    curr2 = curr2->next;
  }
  windowSort(S, 0, S.size(), 2*d, 1, p);
  for (int i = 0; i < S.size(); i++){
    output.push_back(S[i]);
  }
  for (int i = l, j = 0; i < r && j < output.size(); i++, j++) {
    A[i] = output[j];
  }
}

void windowMergeSort(std::vector<int> &A, int l, int r, double p, int d) {
  if (r-l <= 6*d) {
    windowSort(A, l, r, 4*d, d, p);
  }
  else {
    int m = l + (r-l) / 2;
    windowMergeSort(A, l, m, p, d);
    windowMergeSort(A, m, r, p, d);
    windowMerge(A, l, m, r, p, d);
  }
}

std::vector<int> windowOddEvenMerge(std::vector<int> &A, std::vector<int> A1, std::vector<int> A2, double p, double d) {
  if (A1.size() + A2.size() <= 6*d) {
    A1.insert(A1.end(), A2.begin(), A2.end());
    windowSort(A1, 0, int(A1.size()), 4*d, d, p);
    return A1;
  }
  else {
    std::vector<int> A_1o, A_1e, A_2o, A_2e;
    for (int i = 0; i < A1.size(); i++) {
      if (i%2) {
        A_1o.push_back(A1[i]);
      } else {
        A_1e.push_back(A1[i]);
      }
    }
    for (int i = 0; i < A2.size(); i++) {
      if (i%2) {
        A_2o.push_back(A2[i]);
      } else {
        A_2e.push_back(A2[i]);
      }
    }
    std::vector<int> B1 = windowOddEvenMerge(A, A_1e, A_2e, p, d);
    std::vector<int> B2 = windowOddEvenMerge(A, A_1o, A_2o, p, d);
    std::vector<int> B;
    for (int i = 0; i < B1.size(); i++) {
      B.push_back(B1[i]);
      if (i < B2.size()) {
        B.push_back(B2[i]);
      }
    }
    for (int i = 0; i <= B.size()/d; i++) {
      windowSort(B, i*d, std::min<double>(i*d + 6*d, B.size()), 4*d, d, p);
    }
    return B;
  }
}

void windowOddEvenSort(std::vector<int> &A, int l, int r, double p, int d) {
  if (r-l <= 6*d) {
    windowSort(A, l, r, 4*d, d, p);
  }
  else {
    int m = l + (r-l) / 2;
    windowOddEvenSort(A, l, m, p, d);
    windowOddEvenSort(A, m, r, p, d);
    std::vector<int> A1 = std::vector<int>(A.begin()+l, A.begin()+m);
    std::vector<int> A2 = std::vector<int>(A.begin()+m, A.begin()+r);
    std::vector<int> B = windowOddEvenMerge(A, A1, A2, p, d);
    for (int i = l, j = 0; i < r && j < B.size(); i++, j++) {
      A[i] = B[j];
    }
  }
}


std::vector<std::vector<int> > partition(std::vector<int> v){
  std::vector<std::vector<int> > res;
  std::vector<int> shuffled_v = random_sequence_generator(v);
  int sz = sqrt(n);
  int k = log2(n)/2;
  res.push_back(std::vector<int> (v.begin(), v.begin()  + sz));
  int pos = sz;
  int lst;

  for(int i = 1; i <= k; i++){
    lst= std::min(sz + pos, int(shuffled_v.size()));
    res.push_back(std::vector<int> (shuffled_v.begin()+pos, shuffled_v.begin()  + lst));
    pos = pos+sz;
    sz *= 2;
  }
  if(lst<shuffled_v.size()){
    res.push_back(std::vector<int> (shuffled_v.begin()+lst, shuffled_v.begin()  + v.size()));
  }


  return res;
}

std::vector<int> merge(std::vector<int> a, std::vector<int> b, int d, double p){
  std::vector<std::pair<int, int> > rank;
  for(int i = 0; i < b.size(); i++){
    int rank_i = noisyBinarySearch(a, b[i], d, p);
    rank.push_back(std::make_pair(rank_i, b[i]));
  }
  std::sort(rank.begin(), rank.end());
  int i = 0, j = 0;
  std::vector<int> result;
  while(i<a.size() && j<rank.size()){
    if(i<rank[j].first){
      result.push_back(a[i]);
      i++;
    }
    else{
      result.push_back(rank[j].second);
      j++;
    }
  }
  while(i<a.size()){
    result.push_back(a[i]);
    i++;
  }
  while(j<rank.size()){
    result.push_back(rank[j].first);
    j++;
  }
  return result;
}



void riffleSort(std::vector<int> &A, double p) {
  int alpha = 10;
  int gamma = std::max(202*alpha, 909);
  std::vector<std::vector<int> > part = partition(A);
  std::vector<int>  S = part[0];
  windowSort(S, 0, S.size(), std::sqrt(n), 1, p);
  for (int i=0; i<part.size()-1; i++) {
    int d = getMaximumDisloc(S);
    S = merge(S, part[i+1], d, p);
    windowSort(S, 0, S.size(), gamma*9*std::log(n), 1, p);
  }
  A = S;
}



void outputFinalStats(std::string dir_path, long duration) {
  long long inv_count=0, sum_disloc = 0, max_disloc = 0;
  std::unordered_set<long long> seen;
  for(long long i=0;i<output.size();i++){
    if (seen.count(output[i])) {
      std::cout << "duplicate" << std::endl;
    } else {
      seen.insert(output[i]);
    }
    max_disloc = std::max(max_disloc, std::abs(output[i] - i));
    sum_disloc += std::abs(output[i] - i);
    for(int j=i+1;j<output.size();j++) {
      if(output[i]>output[j]) {
        inv_count++;
      }
    }
  }
  for(int i=0;;i++){
    if(exists(dir_path+"/logs_"+std::to_string(i)+".txt"))
      continue;
    else{
      std::ofstream fileout;
      fileout.open(dir_path + "/logs_"+std::to_string(i)+".txt");
      fileout<<"inv_cnt, dis_max, dis_sum, duration, comparisons"<<std::endl;
      fileout<<inv_count<<" "<<max_disloc<<" "<<sum_disloc<<" " << duration << " " <<num_comparisons << std::endl;
      fileout.close();
      break;
    }
  }
}


int main(int argc, char** argv){
  srand(time(NULL));

  fill_hash_table();
  if(argc >= 5){
    n = atoi(argv[3]);
    double p = atof(argv[4]);


    //Creating a directory
    std::string dir_path;
    dir_path = "input_" + std::string(argv[1]) + "_algorithm_" +std::string(argv[2]) + "_n_"+ std::string(argv[3]) +"_p_" +std::string(argv[4]);
    const char * cpath = dir_path.c_str();
    if (mkdir(cpath, 0777) == -1)
      std::cerr << "Error :  " << strerror(errno) << std::endl;

    else
      std::cout << "Directory created" << std::endl;;

    for(int i=0;i<=n-1;i++){
      input.push_back(i);
    }

    if(std::string(argv[1])=="random") {
      random_sequence_generator(input);
    }

    output = input;

    std::chrono::steady_clock::time_point begin, end;
    if(std::string(argv[2])=="insertion"){
      begin = std::chrono::steady_clock::now();
      insertion_sort(output, p);
      end = std::chrono::steady_clock::now();
    }

    if(std::string(argv[2])=="shell") {
      begin = std::chrono::steady_clock::now();
      randomizedShellSort(output, p);
      end = std::chrono::steady_clock::now();
    }

    if(std::string(argv[2])=="shellNo2s3s"){
      begin = std::chrono::steady_clock::now();
      randomizedShellSortWithout2s3s(output, p);
      end = std::chrono::steady_clock::now();
    }

    if(std::string(argv[2])=="shellPratt"){
      begin = std::chrono::steady_clock::now();
      shellSortPratt(output, p);
      end = std::chrono::steady_clock::now();
    }

    if(std::string(argv[2])=="quick"){
      begin = std::chrono::steady_clock::now();
      quickSort(output, 0, n-1, p);
      end = std::chrono::steady_clock::now();
    }

    if(std::string(argv[2])=="annealing"){
      begin = std::chrono::steady_clock::now();
      annealingSort(output, p, 0, 1);
      end = std::chrono::steady_clock::now();
    }

    if(std::string(argv[2])=="riffle"){
      begin = std::chrono::steady_clock::now();
      riffleSort(output, p);
      end = std::chrono::steady_clock::now();
    }

    if(std::string(argv[2])=="windowMerge"){
      begin = std::chrono::steady_clock::now();
      windowMergeSort(output, 0, output.size(), p, std::log2(output.size()));
      end = std::chrono::steady_clock::now();
    }

    if(std::string(argv[2])=="windowOddEven"){
      begin = std::chrono::steady_clock::now();
      windowOddEvenSort(output, 0, output.size(), p, std::log2(output.size()));
      end = std::chrono::steady_clock::now();
    }

    if(std::string(argv[2])=="window"){
      begin = std::chrono::steady_clock::now();
      windowSort(output, 0, output.size(), output.size(), std::log2(output.size()), p);
      end = std::chrono::steady_clock::now();
    }

    auto timeDifference = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();
    outputFinalStats(dir_path, timeDifference);

  }


  return 0;
}

