#ifndef _SIMPLE_COUNTER_
#define _SIMPLE_COUNTER_

#include <iostream>
#include <map>
#include <queue>
#include <algorithm>


using namespace std;

template <class T>
class SimpleCount
{
public:
    SimpleCount(size_t totrack):maxsize(totrack){}
    
    template<typename InputIterator>
    void PushAll(InputIterator first, InputIterator last)
    {
        while (first!=last)
        {
            Push(*first);
            ++first;
        }
    }
    
    void Push(T item)
    {
        auto fentry = mapmajorityset.find(item);
        if( fentry != mapmajorityset.end())
        {
            ++(fentry->second);
        }
        else
        {
            mapmajorityset.insert(make_pair(item,1));
        }
        
        if(mapmajorityset.size() >= maxsize)
        {
            vector< decltype(fentry) > victims;
            
            for (auto mit=mapmajorityset.begin(); mit!=mapmajorityset.end(); ++mit)
            {
                --(mit->second); // first decrement
                
                // map::erase = Iterators, pointers and references referring to elements removed by the function
                // are invalidated. All other iterators, pointers and references keep their validity.
                if(mit->second == 0) { victims.push_back(mit); }
            }
            for(auto qit=victims.begin(); qit!=victims.end(); ++qit)
            {
                mapmajorityset.erase(*qit);
            }
        }
    }
    void PrintAll()
    {
        cout << "We have " << mapmajorityset.size() << " entries " << endl;
        for (auto mit=mapmajorityset.begin(); mit!=mapmajorityset.end(); ++mit)
        {
            cout << mit->first << "(" << mit->second << ") ";
        }
        cout << endl;
    }
    void PrintTop(int currank)
    {
        cout << currank << " has " << mapmajorityset.size() << " entries but only printing top 25" << endl;
        vector< pair<int,T> > majoritysorted;
        for (auto mit=mapmajorityset.begin(); mit!=mapmajorityset.end(); ++mit)
        {
            majoritysorted.push_back(make_pair(mit->second, mit->first));
        }
        sort(majoritysorted.begin(), majoritysorted.end(), std::greater<pair<int,T>>());
        
        int i=0;
        for (auto mit=majoritysorted.begin(); mit!=majoritysorted.end() && i < 25; ++mit, ++i)
        {
            cout << mit->second << "(" << mit->first << ") ";
        }
        cout << endl;
    }
    
    void MergeableSummary(vector<T> & keys, vector<int> & counts)
    {
        for_each(mapmajorityset.begin(),mapmajorityset.end(),
                 [&counts, &keys](pair<const T,int> & p){ counts.push_back(p.second); keys.push_back(p.first); });
    }
    
    bool IsMember(T item)
    {
        return (mapmajorityset.find(item) != mapmajorityset.end());
    }
    
    void CreateIndex()
    {
        for_each(mapmajorityset.begin(),mapmajorityset.end(),
                 [this](pair<const T,int> & p){ internal_.push_back(p); });

    }
    void Clear()
    {
        mapmajorityset.clear();
        vector<pair<T, int>>().swap(internal_);
    }
    
    //! make sure CreateIndex() is called first
    size_t FindIndex(T item) // finds the location of the item in internal_ vector, given the value
    {
        int dummy = 0;
        typename vector<pair<T, int>>::iterator i;
        i = std::lower_bound(internal_.begin(), internal_.end(), make_pair(item,dummy));
        
        size_t ret;
        if (i != internal_.end() && (item == i->first))
            ret = (i - internal_.begin());
        else
            ret = maxsize;
        return ret;
    }
    
    //! make sure CreateIndex() is called first
    T Get(size_t index) // finds the location of the item in internal_ vector, given the value
    {
        return internal_[index].first;
    }
    
    size_t Size()
    {
        return mapmajorityset.size();
    }
    
    void Assign(T * lhskeys, int * lhscounts, size_t lhssize)   // inverse of MergeableSummary
    {
        mapmajorityset.clear(); // delete existing map but retain the maxsize
        for(int i=0; i< lhssize; ++i)
        {
            mapmajorityset.insert(make_pair(lhskeys[i], lhscounts[i]));
        }
    }

    
    void Merge(T * lhskeys, int * lhscounts, size_t lhssize)
    {
        int minweight;
        for (int i = 0; i < lhssize; i++)
        {
            while(lhscounts[i] > 0)
            {
                auto it = mapmajorityset.find(lhskeys[i]);
                if (it != mapmajorityset.end()) // X[z] is the monitored element of a counter c
                {
                    it->second = it->second + lhscounts[i];
                    lhscounts[i] = 0;   // done with this element
                }
                else    // not monitored yet
                {
                    if(mapmajorityset.size() < maxsize)    // there is space for him
                    {
                        mapmajorityset.insert(make_pair(lhskeys[i],lhscounts[i]));
                        lhscounts[i] = 0;   // done with this element
                    }
                    else
                    {
                        // we have to find minimum weight (count) not the minimum key, which is linear if implemented with minelement
                        // however, this is OK because the "else" clause potentially touches every element in the list anyway
                        auto minentry = min_element(mapmajorityset.begin(), mapmajorityset.end(),
                                        [] (typename map<T,int>::value_type & l, typename map<T,int>::value_type & r) { return l.second < r.second; });
                        minweight = minentry->second;
                        
                        int decrementval = std::min(minweight, lhscounts[i]);   // because counts can not go negative
                        vector< typename map<T,int>::iterator > victims;
                        
                        // this process is identical to processing all the "decrementval" entries one after another in a loop
                        for (auto mit=mapmajorityset.begin(); mit!=mapmajorityset.end(); ++mit)
                        {
                            mit->second = mit->second - decrementval; // first decrement
                            if(mit->second == 0) { victims.push_back(mit); }
                        }
                        for(auto qit=victims.begin(); qit!=victims.end(); ++qit)
                            mapmajorityset.erase(*qit);
                        lhscounts[i] -= decrementval;
                    }
                }
            }
		}
	}

    map<T,int> mapmajorityset;
    size_t maxsize;
private:
    vector<pair<T, int>> internal_;
};




#endif
