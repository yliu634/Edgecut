//
// Created by LIU YI on 2021/6/3.
//

#ifndef EDGECUT_LATENCY_ST_PL_H
#define EDGECUT_LATENCY_ST_PL_H

#endif //EDGECUT_LATENCY_ST_PL_H

#pragma once
#include "../common.h"


template<class K, class V, uint IP = 32>
class EdgeCutCommon{
public:
    static const K nn = 100000;
    static const K mm = 4 * nn;
};

void LatencyMatrix(vector<vector<uint16_t>> &la, const bool Seattle = true) {

  uint16_t ns;
  ifstream infile;
  float time;

  if (Seattle) {

    ns = 99;
    vector<uint16_t> lt(ns);
    for (uint16_t i = 0; i < ns; i ++) {
      la.push_back(lt);
    }

    infile.open("../SeattleData_11.txt");
    assert(infile);

    for(uint16_t i = 0; i < ns; i++) {
      for (uint16_t j = 0; j < ns; j++) {
        // while (!infile.eof()) {
        infile >> time;
        time = round(time * 1000); //transfer s to ms;
        la[i][j] = (uint16_t) time;
      }
    }
  } else {
    ns = 490;
    vector<uint16_t> lt(ns);
    for (uint16_t i = 0; i < ns; i ++) {
      la.push_back(lt);
    }

    infile.open("../PlanetLabData_11.txt");
    assert(infile);

    for(uint16_t i = 0; i < ns; i++) {
      for (uint16_t j = 0; j < ns; j++) {
        // while (!infile.eof()) {
        infile >> time;
        time = round(time); //transfer s to ms;
        la[i][j] = time;
      }
    }

  }

  infile.close();
  //std::cout << la.size() << std::endl;

}//end

template<class K, class V, uint IP = 32>
class PreUserList: EdgeCutCommon<K, V, 32>{
public:

    using EdgeCutCommon<K, V, 32>::nn;
    using EdgeCutCommon<K, V, 32>::mm;

    const string str = "../userdataset.txt";

    uint16_t ns;

    PreUserList(){
    }

    explicit PreUserList(const bool Seattle){
      if (Seattle)
        ns = 98;
      else
        ns = 489;
      MobileUserList();
      StaticUserList();
    }

    void MobileUserList(){
      K i(0), k(0); V v;

      ofstream outfile;
      ifstream infile;

      infile.open(str, ios::in);
      assert(infile);

      outfile.open ("../mobileuserlist.txt");
      assert(outfile);

      while(!infile.eof() && i < nn){
        infile >> k; infile >> v;
        outfile << to_string(k) << "\t" << to_string(rand() % ns + 1) << "\n";
        i ++;
      }

      infile.close();
      outfile.close();


    }

    void StaticUserList(){

      K i(0), j(1), k(0); V v;
      //ni = mm / ns;

      ifstream infile;
      ofstream outfile;

      infile.open(str, ios::in);
      assert(infile);

      outfile.open ("../staticuserlist.txt");
      assert(outfile);

      while(!infile.eof() && i < (nn + mm)) {
        infile >> k; infile >> v;
        if (i < nn) {
          i ++;
          continue;
        }

        outfile << to_string(k) << "\t" << to_string(j) << "\n";
        j ++;
        if (j >= ns + 1)
          j = 1;
        i ++;
      }

      infile.close();
      outfile.close();

    }

    void TestKeysGen1(vector<K> &tk, vector<V> &tv, uint32_t lookUpCnt = 1000, const bool Uniform = true) {
      tk.clear();
      tv.clear();

      vector <K> keys(nn);
      vector <V> values(nn);
      K i(0), idx(0);

      ifstream infile;
      infile.open("../mobileuserlist.txt", ios::in);
      assert(infile);

      while (!infile.eof() && i < nn) {
        infile >> keys[i];
        infile >> values[i];
        i ++;
      }

      i = 0;
      infile.close();

      if (Uniform) {
        srand(0);
        for (i = 0; i < lookUpCnt; i++) {
          idx = rand() % nn;
          tk.push_back(keys[idx]);
          tv.push_back(values[idx]);
        }
      } else {
        for (i = 0; i < lookUpCnt; i++) {
          InputBase::distribution = zipfian;
          InputBase::bound = nn;
          InputBase::setSeed(time(NULL));
          idx = InputBase::rand();
          tk.push_back(keys[idx]);
          tv.push_back(values[idx]);
        }
      }

    }

    void TestKeysGen2(vector<K> &tk, vector<V> &tv, uint32_t lookUpCnt = 1000, const bool Uniform = true, const bool Static = true){
      tk.clear();
      tv.clear();

      vector <K> keys(nn);
      vector <V> values(nn);

      vector <K> statickeys(mm);
      vector <V> staticvalues(mm);

      K i(0), idx(0);

      uint32_t nnmm = (uint32_t)(round( lookUpCnt / ((mm/nn) + 1)));

      ifstream infile;
      infile.open("../mobileuserlist.txt", ios::in);
      assert(infile);

      while (!infile.eof() && i < nn) {
        infile >> keys[i];
        infile >> values[i];
        i ++;
      }

      i = 0;
      infile.close();

      // again
      {
        infile.open("../staticuserlist.txt", ios::in);
        assert(infile);

        while (!infile.eof() && i < nn) {
          infile >> statickeys[i];
          infile >> staticvalues[i];
          i++;
        }

        i = 0;
        infile.close();
      }


      if (Uniform) {
        srand(0);
        for (i = 0; i < nnmm ; i++) {
          idx = rand() % nn;
          tk.push_back(keys[idx]);
          tv.push_back(values[idx]);
        }

        for (i = nnmm; i < lookUpCnt; i++) {
          idx = rand() % mm;
          tk.push_back(statickeys[idx]);
          tv.push_back(staticvalues[idx]);
        }

      } else {

        InputBase::distribution = zipfian;
        InputBase::bound = nn;
        InputBase::setSeed(time(NULL));

        for (i = 0; i < nnmm; i++) {

          idx = InputBase::rand();
          tk.push_back(keys[idx]);
          tv.push_back(values[idx]);
        }

        for (i = nnmm; i < lookUpCnt; i++) {
          idx = rand() % mm;
          tk.push_back(statickeys[idx]);
          tv.push_back(staticvalues[idx]);
        }

      }

    }

    ~PreUserList(){

    }

};

