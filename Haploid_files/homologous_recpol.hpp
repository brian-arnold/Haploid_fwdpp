#ifndef HOMOLOGOUS_REC_HPP__
#define HOMOLOGOUS_REC_HPP__

#include <vector>
#include <algorithm>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <unordered_set>


struct homologous_rec_SLOW
/*!
  \brief Simple model of homologous recombination.
  Generates a Poisson-distributed number of recombination breakpoints with
  mean recrate that are uniformly-distributed between minpos and maxpos,
 then for each breakpoint, a second one is added via a geometric dist
 to simulate short gene conversion tracts
 */
{
    const double recrate, minpos, maxpos, meantrlen, fragsize, geop ;
    mutable std::vector<std::pair<double,double>> rawpos;
    // mutable to allow it to be modified by const func below
    
    explicit homologous_rec_SLOW(const double recrate_, const double minpos_,
                           const double maxpos_, const double meantrlen_,
                            const double fragsize_)
        : recrate{ recrate_ }, minpos{ minpos_ }, maxpos{ maxpos_ },
            meantrlen{ meantrlen_ }, fragsize{ fragsize_ },
            geop{ 1.0/meantrlen_ }
    /*!
      \param recrate_ The recombination rate (per diploid_, per region)
      \param minpos_ The minimum recombination position allowed
      \param maxpos_ The maximum recombination position allowed
      the gametes & mutations involve in an x-over.
     */
    {
    }
    
    homologous_rec_SLOW(const homologous_rec_SLOW &) = default;
    homologous_rec_SLOW(homologous_rec_SLOW &&) = default;
    
    std::vector<double>
    operator()(const gsl_rng * r) const
    {
        rawpos.clear() ;
        
        
        
        
        
        // calculate number of recomb. events
        unsigned nbreaks
            = (recrate > 0) ? gsl_ran_poisson(r, recrate) : 0u;
        if (!nbreaks){
            //std::cout << "no recs!\n" ;
            return {};
        }
        
        //std::cout << "recs:\t" << nbreaks << "\n" ;
      /*
        std::cout << "geop: " << geop << "\n";
        for (unsigned i = 0; i < nbreaks; ++i){
            std::cout << gsl_ran_geometric(r, geop)/fragsize << "\n" ;
        }
        */
        rawpos.reserve(nbreaks) ;
        std::vector<double> pos ;
        pos.reserve(nbreaks + 1);
        std::vector<int> toDelete ;
        // generate breakpoints for rec events
        for (unsigned i = 0; i < nbreaks; ++i){
            // returns a random variate from the flat (uniform) distribution
            double recpt = gsl_ran_flat(r, minpos, maxpos) ;
            if (gsl_rng_uniform(r) < 0.5){
                // rec proceeds to right, no breakpoints beyond fragsize
                // geometric distribution p(1-p)^(k-1), which has
                // mean 1/p
                rawpos.push_back(std::make_pair(recpt,
                                                std::min((recpt+(gsl_ran_geometric(r, geop)/fragsize)),
                                                maxpos))) ;
                
                //pos.emplace_back( recpt ) ;
                //pos.emplace_back( std::min((recpt+(gsl_ran_geometric(r, geop)/fragsize)),
                //                           maxpos) ) ;
            }else{
                // rec proceeds to left, no negative breakpoints
                rawpos.push_back(std::make_pair(std::max((recpt-(gsl_ran_geometric(r, geop)/fragsize)),
                                                        minpos), recpt)) ;
                
                //pos.emplace_back( recpt ) ;
                //pos.emplace_back( std::max((recpt-(gsl_ran_geometric(r, geop)/fragsize)),
                //                           minpos) ) ;
            }
        }
        //std::sort(rawpos.begin(), rawpos.end()) ;
        // test if 2+ rec events have overlapping boundaries, if so merge them
        /*
         *** Types of overlapping recombination events ***
         chromo: -------------------------------------------------------
         rec1:                   --------
         rec2:                   ------------               * ends equal, likely at chromo boundary 0
         rec3:                        ------------
         rec4:                           ------------       * ends equal
         rec5:               ------------                   * ends equal, likely at chromo boundary 1
         rec6:          ------------
         rec7:       ------------                           * ends equal
         rec8:                      ---
         rec9:                   --------                   * both ends equal
         rec10:           ----------------------
         rec11:   ----                                      * no conflict
         rec12:                                       ----  * no conflict
        comparing rec event A with event B (e.g. A=rec1, B=rec2-12
             if B.min > A.max || B.max < A.min
                #majority of cases, no conflicts
                #rec1/rec11, rec1/rec12, respectively
             else
                 if B.min <= A.max && B.min >= A.min        //B.min within interval of A
                    #rec1 conflicts with rec2/3/4
                    #set A.max=B.max, delete B
         
                 if B.max >= A.min && B.max <= A.max        //B.max within interval of A
                    #rec1 conflicts with rec5/6/7
                    #set A.min=B.min, delete B
         
                 if B.min >= A.min && B.max <= A.max        //B subset or identical to A
                    #rec1 conflict with rec8/9
                    #delete B
         
                 if B.min < A.min && B.max > A.max          //A subset of B
                    #rec1 conflict with rec10
                    #delete A
         */
        /*
        std::cout << "#############b4 conflict resolution\t" <<  rawpos.size() << "\n";
        for(unsigned i = 0; i<nbreaks; ++i){
            std::cout << rawpos[i].first << "\t" << rawpos[i].second << "\n" ;
        }
         */
        if(nbreaks>1){
            for (unsigned a = 0; a < nbreaks-1; a++){
                for (unsigned b = a+1; b < nbreaks; b++){
                    if(rawpos[b].first > rawpos[a].second || rawpos[b].second < rawpos[a].first){
                        // do nothing
                    }else{
                        // the order of these statements matter!!
                        // first check if one rec event is subset of other
                        // if not, then rec events staggered, with 1 brkpoint in interval
                        if(rawpos[b].first >= rawpos[a].first && rawpos[b].second <= rawpos[a].second){
                            //rec1 conflict with rec8/9
                            //delete B
                            toDelete.push_back(b) ;
                        }else if(rawpos[b].first < rawpos[a].first && rawpos[b].second > rawpos[a].second){
                            //rec1 conflict with rec10
                            //delete A
                            toDelete.push_back(a) ;
                        }else if(rawpos[b].first <= rawpos[a].second && rawpos[b].first >= rawpos[a].first){
                            //rec1 conflicts with rec2/3/4
                            //set A.max=B.max, delete B
                            rawpos[a].second = rawpos[b].second ;
                            toDelete.push_back(b) ;
                        }else if(rawpos[b].second >= rawpos[a].first && rawpos[b].second <= rawpos[a].second){
                            //rec1 conflicts with rec5/6/7
                            //set A.min=B.min, delete B
                            rawpos[a].first = rawpos[b].first ;
                            toDelete.push_back(b) ;
                        }
                    }
                }
            }
            // sort in reverse order to delete from vector going in reverse
            std::sort( toDelete.rbegin(), toDelete.rend() );
            // unique returns iterator to new last element, erase everything after
            toDelete.erase(std::unique(toDelete.begin(), toDelete.end()),
                           toDelete.end()) ;
            /*
            std::cout << "#############Right before deleting:\n" ;
            for(unsigned i = 0; i<nbreaks; ++i){
                std::cout << rawpos[i].first << "\t" << rawpos[i].second << "\n" ;
            }
            std::cout << "\n" ;
            */
            for(unsigned d = 0; d<toDelete.size() ; ++d){
                rawpos[toDelete[d]] = rawpos.back() ;
                rawpos.pop_back() ;
            }
           
            //std::cout << "#############aft conflict resolution\t" <<  rawpos.size() << "\n"; ;
            for(unsigned i = 0; i<rawpos.size(); ++i){
                //std::cout << rawpos[i].first << "\t" << rawpos[i].second << "\n" ;
            }
        }
        
        for(unsigned i = 0; i<rawpos.size(); ++i){
            //std::cout << rawpos[i].first << "\t" << rawpos[i].second << "\n" ;
            pos.push_back(rawpos[i].first) ;
            pos.push_back(rawpos[i].second) ;
        }
        std::sort(pos.begin(), pos.end());
        // Note: this is required for all vectors of breakpoints!
        pos.push_back(std::numeric_limits<double>::max());
        /*
        for(unsigned i = 0; i<pos.size(); ++i){
            std::cout << pos[i] << "\t" ;
        }
        std::cout << "\n" ;
        */
        return pos;
    }
};


struct homologous_rec
/*!
 \brief Simple model of homologous recombination.
 Generates a Poisson-distributed number of recombination breakpoints with
 mean recrate that are uniformly-distributed between minpos and maxpos,
 then for each breakpoint, a second one is added via a geometric dist
 to simulate short gene conversion tracts
 */
{
    const double recrate, minpos, maxpos, meantrlen, fragsize, geop ;
    // pos is a vector of break points, ocurring in pairs
    // e.g. pos[0] and pos[1] are the begin and end of rec1
    // mutable to allow it to be modified by const func below
    mutable std::vector<double> pos;
    
    explicit homologous_rec(const double recrate_, const double minpos_,
                            const double maxpos_, const double meantrlen_,
                            const double fragsize_)
    : recrate{ recrate_ }, minpos{ minpos_ }, maxpos{ maxpos_ },
    meantrlen{ meantrlen_ }, fragsize{ fragsize_ },
    geop{ 1.0/meantrlen_ }
    /*!
     \param recrate_ The recombination rate (per diploid_, per region)
     \param minpos_ The minimum recombination position allowed
     \param maxpos_ The maximum recombination position allowed
     the gametes & mutations involve in an x-over.
     */
    {
    }
    
    homologous_rec(const homologous_rec &) = default;
    homologous_rec(homologous_rec &&) = default;
    
    std::vector<double>
    operator()(const gsl_rng * r) const
    {
        pos.clear() ;

        // calculate number of recomb. events
        unsigned nbreaks
            = (recrate > 0) ? gsl_ran_poisson(r, recrate) : 0u;
        if (!nbreaks){
            //std::cout << "no recs!\n" ;
            return {};
        }
        
        //std::cout << "recs:\t" << nbreaks << "\n" ;
        /*
         std::cout << "geop: " << geop << "\n";
         for (unsigned i = 0; i < nbreaks; ++i){
         std::cout << gsl_ran_geometric(r, geop)/fragsize << "\n" ;
         }
         */
        std::vector<double> pos ;
        pos.reserve(2*nbreaks + 1);
        std::vector<int> toDelete ;
        // generate breakpoints for rec events
        for (unsigned i = 0; i < nbreaks; ++i){
            // returns a random variate from the flat (uniform) distribution
            double recpt = gsl_ran_flat(r, minpos, maxpos) ;
            if (gsl_rng_uniform(r) < 0.5){
                // rec proceeds to right, no breakpoints beyond fragsize
                // geometric distribution p(1-p)^(k-1), which has
                // mean 1/p
                pos.push_back(recpt) ;
                pos.push_back(std::min((recpt+(gsl_ran_geometric(r, geop)/fragsize)),
                                        maxpos)) ;
            }else{
                // rec proceeds to left, no negative breakpoints
                pos.push_back(std::max((recpt-(gsl_ran_geometric(r, geop)/fragsize)),
                                       minpos)) ;
                pos.push_back(recpt) ;
            }
        }
        //std::sort(rawpos.begin(), rawpos.end()) ;
        // test if 2+ rec events have overlapping boundaries, if so merge them
        /*
         *** Types of overlapping recombination events ***
         chromo: -------------------------------------------------------
         rec1:                   --------
         rec2:                   ------------               * ends equal, likely at chromo boundary 0
         rec3:                        ------------
         rec4:                           ------------       * ends equal
         rec5:               ------------                   * ends equal, likely at chromo boundary 1
         rec6:          ------------
         rec7:       ------------                           * ends equal
         rec8:                      ---
         rec9:                   --------                   * both ends equal
         rec10:           ----------------------
         rec11:   ----                                      * no conflict
         rec12:                                       ----  * no conflict
         comparing rec event A with event B (e.g. A=rec1, B=rec2-12
         if B.min > A.max || B.max < A.min
         #majority of cases, no conflicts
         #rec1/rec11, rec1/rec12, respectively
         else
         if B.min <= A.max && B.min >= A.min        //B.min within interval of A
         #rec1 conflicts with rec2/3/4
         #set A.max=B.max, delete B
         
         if B.max >= A.min && B.max <= A.max        //B.max within interval of A
         #rec1 conflicts with rec5/6/7
         #set A.min=B.min, delete B
         
         if B.min >= A.min && B.max <= A.max        //B subset or identical to A
         #rec1 conflict with rec8/9
         #delete B
         
         if B.min < A.min && B.max > A.max          //A subset of B
         #rec1 conflict with rec10
         #delete A
         */
        /*
         std::cout << "#############b4 conflict resolution\t" <<  pos.size() << "\n";
         for(unsigned i = 0; i < 2*nbreaks-1; i+=2){
             std::cout << pos[i] << "\t" << pos[i+1] << "\n" ;
         }
        */
        if(nbreaks>1){
            for (unsigned a = 0; a < 2*nbreaks-2; a+=2){
                for (unsigned b = a+2; b < 2*nbreaks-1; b+=2){
                    if(pos[b] > pos[a+1] || pos[b+1] < pos[a]){
                        // do nothing
                    }else{
                        // the order of these statements matter!!
                        // first check if one rec event is subset of other
                        // if not, then rec events staggered, with 1 brkpoint in interval
                        // not that a and b are starting points of recs, a+1 and b+1 are endpoints
                        if(pos[b] >= pos[a] && pos[b+1] <= pos[a+1]){
                            //rec1 conflict with rec8/9
                            //delete B
                            toDelete.push_back(b) ;
                        }else if(pos[b] < pos[a] && pos[b+1] > pos[a+1]){
                            //rec1 conflict with rec10
                            //delete A
                            toDelete.push_back(a) ;
                        }else if(pos[b] <= pos[a+1] && pos[b] >= pos[a]){
                            //rec1 conflicts with rec2/3/4
                            //set A.max=B.max, delete B
                            pos[a+1] = pos[b+1] ;
                            toDelete.push_back(b) ;
                        }else if(pos[b+1] >= pos[a] && pos[b+1] <= pos[a+1]){
                            //rec1 conflicts with rec5/6/7
                            //set A.min=B.min, delete B
                            pos[a] = pos[b] ;
                            toDelete.push_back(b) ;
                        }
                    }
                }
            }
            // sort in reverse order to delete from vector going in reverse
            std::sort( toDelete.rbegin(), toDelete.rend() );
            // unique returns iterator to new last element, erase everything after
            toDelete.erase(std::unique(toDelete.begin(), toDelete.end()),
                           toDelete.end()) ;
            /*
             std::cout << "#############Right before deleting:\n" ;
            for(unsigned i = 0; i < 2*nbreaks-1; i+=2){
                std::cout << pos[i] << "\t" << pos[i+1] << "\n" ;
            }
            std::cout << "\n" ;
             */
            for(unsigned d = 0; d<toDelete.size() ; ++d){
                pos[toDelete[d]+1] = pos.back() ;// delete endpoint
                pos.pop_back() ;
                pos[toDelete[d]] = pos.back() ;// delete starting point
                pos.pop_back() ;
            }
            /*
            //std::cout << "#############aft conflict resolution\t" <<  rawpos.size() << "\n"; ;
            std::cout << "final positions:\t" << pos.size() << "\n" ;
           for(unsigned i = 0; i<pos.size()-1; i+=2){
                std::cout << pos[i] << "\t" << pos[i+1] << "\n" ;
            }
             */
        }
        
        std::sort(pos.begin(), pos.end());
        // Note: this is required for all vectors of breakpoints!
        pos.push_back(std::numeric_limits<double>::max());
        /*
         for(unsigned i = 0; i<pos.size(); ++i){
         std::cout << pos[i] << "\t" ;
         }
         std::cout << "\n" ;
        */
        return pos;
    }
};



#endif
