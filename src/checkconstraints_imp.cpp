#include <cmath>
#include "checkconstraints.h"

//1 = head/householder, 2 = spouse, 3 = child, 4 = child-in-law, 5 = parent, 6 = parent-in- law, 7 = sibling, 8 = sibling-in-law,
//9 = grandchild,
//10 = other relatives,
//11 = partner, friend, visitor,
//12 = other non-relatives
inline bool IsHead(double relate, double age) {
    return (relate == HEAD && age >=16);
}

inline int GetHead(double *record, int hhsize) {
    for (int i = 1; i <= hhsize; i++) {
        if (record[2*hhsize+i]==HEAD) {
            return i;
        }
    }
    return -1;
}

inline bool MoreThanOneHead(double *record, int hhsize) {
    int nhead = 0;
    for (int i = 1; i <= hhsize; i++) {
        if (record[2*hhsize+i]==HEAD) {
            nhead++;
        }
    }
    return (nhead >1);
}

inline int GetValidSpouse(double *record, int hhsize) {
    int nspouse = 0;
    int spouse = -1;
    for (int i = 1; i <= hhsize; i++) {
        if (record[2*hhsize+i]==SPOUSE) {
            nspouse++;
            spouse = i;
        }
    }
    if (nspouse > 1) {return 0;} //too many spouse
    if (nspouse == 0) { return -1;} //no spouse
    if (record[hhsize+spouse]<16) {return 0;} //spouse is under-age
    return spouse;
}

inline bool IsValidCouple(double *record, int hh_size, int spouse, int head) {
    if (spouse ==0) { //bad spouse or too many spouses
        return false;
    } else { //valid spouse or no spouse
        if (spouse>0) {//the only spouse, so check sex, and age difference
            if (record[head] == record[spouse]) {return false;}
            if (std::abs(record[hh_size + head] - record[hh_size + spouse]) > 49) {return false;}
        }
    }
    return true;
}

//return -1 if no child
//return the record index of the oldest child otherwise
inline int GetOldestChild(double *record, int hhsize) {
    double age = -1;
    int child = -1;  //no childen
    for (int i = 1; i <= hhsize; i++) {
        if (record[2*hhsize+i]==CHILD) {
            if (record[hhsize+i] > age) {
                age = record[hhsize+i];
                child = i;
            }
        }
    }
    return child;
}

inline bool IsValidChild(double *record, int hh_size, int child, int head) {
    if (child>0) {//get a child, check age difference
        if (record[hh_size + head] - record[hh_size + child] <7) {return false;}
    }
    return true;
}


//return -1 if no child
//return the record index of the oldest child otherwise
//inline int GetOldestChildInLaw(double *record, int hhsize) {
//    double age = -1;
//    int child = -1;  //no childen
//    for (int i = 1; i <= hhsize; i++) {
//        if (record[2*hhsize+i]==CHILDINLAW) {
//            if (record[hhsize+i] > age) {
//                age = record[hhsize+i];
//                child = i;
//            }
//        }
//    }
//    return child;
//}

//inline bool IsValidChildInLaw(double *record, int hh_size, int child, int head) {
//    if (child>0) {//get a child, check age difference
//        if (record[hh_size + head] - record[hh_size + child] <10) {return false;}
//    }
//    return true;
//}

//return -1 if no parent
//return the record index of the youngest parent otherwise
inline int GetYoungestParent(double *record, int hhsize) {
    double age = 1000;
    int parent = -1;  //no childen
    for (int i = 1; i <= hhsize; i++) {
        if (record[2*hhsize+i]==PARENT) {
            if (record[hhsize+i] < age) {
                age = record[hhsize+i];
                parent = i;
            }
        }
    }
    return parent;
}

inline bool IsValidParent(double *record, int hh_size, int parent, int head) {
    if (parent>0) {//get a child, check age difference
        if (record[hh_size + parent] -record[hh_size + head] <10) {return false;}
    }
    return true;
}

inline int GetYoungestParentInLaw(double *record, int hhsize) {
    double age = 1000;
    int parent = -1;  //no childen
    for (int i = 1; i <= hhsize; i++) {
        if (record[2*hhsize+i]==PARENTINLAW) {
            if (record[hhsize+i] < age) {
                age = record[hhsize+i];
                parent = i;
            }
        }
    }
    return parent;
}

inline bool IsValidParentInLaw(double *record, int hh_size, int parent, int head) {
    if (parent>0) {//get a child, check age difference
        if (record[hh_size + parent] -record[hh_size + head] <4) {return false;}
    }
    return true;
}

inline bool IsValidSiblingOrSiblingInLaw(double *record, int hhsize, int head) {
    for (int i = 1; i <= hhsize; i++) {
        if (record[2*hhsize+i]== SIBLING || record[2*hhsize+i] == SIBLINGINLAW) {
            if (std::abs(record[hhsize + i] - record[hhsize + head]) >37) {return false;}
        }
    }
    return true;
}

inline bool IsValidGrandChild(double *record, int hhsize, int spouse, int head) {
    for (int i = 1; i <= hhsize; i++) {
        if (record[2*hhsize+i]== GRANDCHILD) {
            if (record[hhsize + head] < 33) {return false;} //too young to be grand parent for the HEAD
            if (spouse > 0) { //make sure the spouse(if any) is not too young
                if (record[hhsize + spouse] < 17) {return false;}
            }
            if (record[hhsize + head] - record[hhsize + i] <26 ) {return false;}
        }
    }
    return true;
}


int isValid(double *datah, int hh_size) {

    int head = GetHead(datah,hh_size);
    if (head <=0) {return 0;}

    if (!IsHead(datah[2 * hh_size + head], datah[hh_size + head])) {return 0;}
    if (MoreThanOneHead(datah,hh_size)) {return 0;}

    int spouse = GetValidSpouse(datah,hh_size);
    if (!IsValidCouple(datah,hh_size,spouse, head)) {return 0;}

    int oldestChild = GetOldestChild(datah,hh_size);
    if (!IsValidChild(datah,hh_size,oldestChild,head)) {return 0;}

    //int oldestChildInLaw = GetOldestChildInLaw(datah,hh_size);
    //if (!IsValidChildInLaw(datah,hh_size,oldestChildInLaw,head)) {return 0;}

    int youngestParent = GetYoungestParent(datah,hh_size);
    if (!IsValidParent(datah,hh_size,youngestParent,head)) {return 0;}

    int youngestParentInLaw = GetYoungestParentInLaw(datah,hh_size);
    if (!IsValidParentInLaw(datah,hh_size,youngestParentInLaw,head)) {return 0;}

    if (!IsValidSiblingOrSiblingInLaw(datah,hh_size,head)) {return 0;}

    if (!IsValidGrandChild(datah,hh_size,spouse,head)) {return 0;}

    return 1;

}

int checkconstraints_imp(double *data, double *isPossible,int hh_size, int DIM, int nHouseholds) {

    int totalpossible = 0;
    double *datah = new double[hh_size * 3 + 1];
    //column 3, 6, 7 = sex, age and relte
    int column[COL]; column[0] = 3; column[1] = 6; column[2] = 7;

	for (int m = 1; m <= nHouseholds; m++){
        for (int j = 1; j <= hh_size; j++) {
            for (int k = 0; k < COL; k++) {
                datah[k * hh_size + j] = data[((j-1) * DIM + column[k] -1) * nHouseholds + (m-1)];
            }
        }
		isPossible[m-1] = isValid(datah, hh_size);
        totalpossible+= (int)isPossible[m-1];
	}

	delete [] datah;
    return totalpossible;
}

int checkconstraints_imp_HHhead_at_group_level(double *data, double *isPossible,int hh_size, int DIM, int nHouseholds) {
  int realsize = hh_size + 1;
  int totalpossible = 0;
  double *datah = new double[realsize * 3 + 1];
  //column 2, 5, 6 = sex, age and relte //zero-based here
  int column[COL]; column[0] = 0; column[1] = 3; column[2] = 4;

  for (int m = 1; m <= nHouseholds; m++){
    //for each household member
    for (int j = 1; j < realsize; j++) {
      for (int k = 0; k < COL; k++) {
        datah[k * realsize + j] = data[((j-1) * DIM + column[k]+2) * nHouseholds + (m-1)];
        if (k+1 == COL) { //relate column
          datah[k * realsize + j] = datah[k * realsize + j] + 1; //addjust by adding 1
        }
      }
    }
    datah[realsize] = data[(column[0]+8) * nHouseholds + (m-1)];
    datah[2 * realsize] = data[(column[1]+8) * nHouseholds + (m-1)];
    datah[3 * realsize] = 1;

    isPossible[m-1] = isValid(datah, realsize);
    totalpossible+= (int)isPossible[m-1];
  }

  delete [] datah;
  return totalpossible;
}
