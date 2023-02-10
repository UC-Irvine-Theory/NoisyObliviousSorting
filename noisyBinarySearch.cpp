#include <iostream>
#include <vector> 
#include <cmath>
#include <random>
#include <algorithm>
#include <memory>
#include <limits>
#include "noisyBinarySearch.hpp"
struct node
{
  int intervalStart;
  int intervalEnd;
  node *left;
  node *right;
  node *parent;
  std::shared_ptr<int> L;
  std::shared_ptr<int> R;
  int leafIndex;
};

class noisyTree
{
public:
  noisyTree(): root(nullptr) {constructNoisyTree();}
  noisyTree(int j, int c, int d, int h, int eta, int k, double p): j(j), c(c), d(d), h(h), eta(eta), k(k), p(p), root(nullptr) {constructNoisyTree();}
  ~noisyTree();
  void constructNoisyTree();
  std::pair<bool, int> walkOnTree(int x, std::vector<int> &S);
  bool test(int x, node *v, std::vector<int> &S);

private:
  void populateNoisyTree(node *node);
  void destroyTree(node *node);
  node *root;
  int j;
  int h;
  int c;
  int d;
  int eta;
  int k;
  double p;
};

noisyTree::~noisyTree()
{
  destroyTree(root);
}

void noisyTree::destroyTree(node *node)
{
  if (node)
  {
    destroyTree(node->left);
    destroyTree(node->right);
    delete node;
  }
}

void noisyTree::populateNoisyTree(node* node) {
  if (node) {
    populateNoisyTree(node->left);
    populateNoisyTree(node->right);
    if (!node->left and !node->right) {
      int groupId = 2*node->leafIndex+j;
      node->intervalStart = k*c*groupId*d ;
      node->intervalEnd = node->intervalStart + k*c*d - 1;
      node->L = std::shared_ptr<int>(new int(node->intervalStart-d-1));
      node->R = std::shared_ptr<int>(new int(node->intervalEnd+d));
    }
    else if (node->left and !node->right) {
      node->intervalStart = node->left->intervalStart;
      node->intervalEnd = node->left->intervalEnd;
      node->L = node->left->L;
      node->R = node->left->R;
    }
    else if (node->left and node->right) {
      node->intervalStart = std::min(node->left->intervalStart, node->right->intervalStart);
      node->intervalEnd = std::max(node->left->intervalEnd, node->right->intervalEnd);
      node->L = node->left->L;
      node->R = node->right->R;
    }
  }
}

// constructs binary tree with first h+1 levels complete, then with 2^h paths of length eta attached to the leaves, for a tree with total height of h+eta
void noisyTree::constructNoisyTree()
{
  // complete binary tree with h+1 levels
  root = new node{-1, -1, nullptr, nullptr, nullptr, nullptr, nullptr, -1};
  std::vector<node*> currentLevel = {root};
  for (int i = 0; i < h; i++) {
    std::vector<node*> newLevel;
    for (node *currNode: currentLevel) {
      currNode -> left = new node{-1, -1, nullptr, nullptr, currNode, nullptr, nullptr, -1};
      currNode -> right = new node{-1, -1, nullptr, nullptr, currNode, nullptr, nullptr, -1};
      newLevel.push_back(currNode->left);
      newLevel.push_back(currNode->right);
    }
    currentLevel = newLevel;
  }
  
  // add paths of length eta to each leaf of the complete tree
  for (int i = 0; i < eta; i++) {
    std::vector<node*> newLevel;
    for (node *currNode: currentLevel) {
      currNode -> left = new node{-1, -1, nullptr, nullptr, currNode, nullptr, nullptr, -1};
      newLevel.push_back(currNode->left);
    }
    currentLevel = newLevel;
  }
  
  // index the leaves of tree
  for (int i = 0; i<currentLevel.size();i++) {
    currentLevel[i]->leafIndex = i;
  }
  
  // populate intervalStart, intervalEnd, L, and R for each node
  populateNoisyTree(root);
}

// first element of pair is true if walk successful, second element is the computed rank
std::pair<bool, int> noisyTree::walkOnTree(int x, std::vector<int> &S) {
  int iterationsRemaining = 7*std::floor(std::log2(S.size()));
  node *curr = root;
  while (iterationsRemaining--) {
    // if curr is a leaf node, need to return any position in the interval of curr. 
    // we just return the start of the interval
    if (curr->leafIndex != -1) {
      return std::make_pair(true, curr->intervalStart);
    }
    
    // test x with both left and right
    bool testL = test(x, curr->left, S), testR = test(x, curr->right, S);
    
    // if only one of them is true, walk to that node
    // if both are false, walk to parent if it exists
    // otherwise, stay at curr
    if (testL and !testR) {
      curr = curr->left;
    } else if (!testL and testR) {
      curr = curr->right;
    } else if (!testL and !testR and curr->parent) {
      curr = curr->parent;
    }
  }
  // if no iterations remaining, walk was unsuccessful. return an arbitrary position of 0 as the result
  return std::make_pair(false, 0);
}

// test operation described in [GLLP'18]
bool noisyTree::test(int x, node *v, std::vector<int> &S) {
  if (!v) {
    return false;
  }
  int failingL = 0, failingR = 0;
  int SL, SR;
  for (int i = 0; i<k; i++) {
    // check if L(v) or R(v) are out of bounds, and if so set SL, SR to -inf,+inf respectively
    if (*(v->L) < 0) {
      SL = std::numeric_limits<int>::min();
    } else {
      SL = S[*(v->L)];
    }

    if (*(v->R) > S.size()-1) {
      SR = std::numeric_limits<int>::max();
    } else {
      SR = S[*(v->R)];
    }

    // decrement L(v) and increment R(v)
    (*(v->L))--;
    (*(v->R))++;
    // test succeeds if both x>SL and x<SR holds, and fails otherwise. If either SL or SR is -inf or +inf respectively, that case always passes
    // since we are assuming that comparisons with +-inf will not be subject to errors
    int idx_SL = std::distance(S.begin(), std::find(S.begin(), S.end(), SL));
    int idx_SR = std::distance(S.begin(), std::find(S.begin(), S.end(), SR));
    int idx_x = std::distance(S.begin(), std::find(S.begin(), S.end(), x));
    failingL += (SL != std::numeric_limits<int>::min() and !is_inversion_with_noise(S, idx_SL, idx_x, p));
    failingR += (SR != std::numeric_limits<int>::max() and !is_inversion_with_noise(S, idx_x, idx_SR, p));
  }

  if (failingL > k/2 or  failingR > k/2) {
    return false;
  } else {
    return true;
  }
}


// compute approximate rank of x in S, given an error probability p
int noisyBinarySearch(std::vector<int> &S, int x, int d, double p) {
  int k = 1;
  int c = 5;
  if (d == 0 or (S.size()+1)/(2*k*c*d) == 0) {
    return 0;
  }
  int h = std::log2((S.size()+1)/(2*k*c*d));
  int eta = 2*std::ceil(std::log2(S.size()));
  noisyTree T0(0, c, d, h, eta, k, p);
  noisyTree T1(1, c, d, h, eta, k, p);

  // bool represents whether walk was successful, and int is the computed rank.
  // if the walk was unsuccessful, an arbitrary rank of 0 is returned
  std::pair<bool, int> r = T0.walkOnTree(x, S);
  if (r.first) {
    return r.second;
  } else {
    r = T1.walkOnTree(x, S);
    return r.second;
  }
}
