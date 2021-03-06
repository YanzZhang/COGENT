#ifndef  _CLSINTERFACE_H_
#define  _CLSINTERFACE_H_

#include <string>

#include "REAL.H"
#include "ParmParse.H"
#include "KineticSpecies.H"
#include "GlobalDOF.H"
#include "MultiPhysicsPreconditioner.H"
#include "NamespaceHeader.H"


class GKVector;
class GKOps;
/**
 * Collision interface class.
 *
 * Abstract base class the provides an interface on which algorithm classes
 * can rely to set collision on species KineticSpecies.
*/

class CLSInterface
{
   public:

      /// Destructor.
      /**
       */
      virtual ~CLSInterface() {;}

      /// Evaluate collison RHS.
      /**
       *  Pure virtual function that will be resposible for evaluating
       *  collision RHS on the gven species.
       *
       *  @param[in,out] rhs - collision contribution to be added to rhs.
       *  @param[in] soln - kinetic species.
       *  @param[in] species - species to which collision model is applied.
       *  @param[in] time the time at which the collison RHS is to be evaluated
       */
      virtual void evalClsRHS( KineticSpeciesPtrVect&       rhs,
                               const KineticSpeciesPtrVect& soln,
                               const int                    species_a,
                               const int                    species_b,
                               const Real                   time ) = 0;

    
      virtual void evalClsRHSImplicit( KineticSpeciesPtrVect&       rhs,
                                       const KineticSpeciesPtrVect& soln,
                                       const int                    species_a,
                                       const int                    species_b,
                                       const Real                   time ) = 0;

      virtual void evalClsRHSExplicit( KineticSpeciesPtrVect&       rhs,
                                       const KineticSpeciesPtrVect& soln,
                                       const int                    species_a,
                                       const int                    species_b,
                                       const Real                   time ) = 0;


      virtual inline bool isLinear() { return false; }

      virtual Real computeDt(const KineticSpeciesPtrVect& soln, const int i, const int j) { return DBL_MAX; }
      virtual Real TimeScale(const KineticSpeciesPtrVect& soln, const int i, const int j) { return -1; }

      virtual int precondMatrixBands() { return(0); }

      virtual void assemblePrecondMatrix( void*,
                                          const KineticSpecies&,
                                          const GlobalDOFKineticSpecies&,
                                          const Real ) { }

      virtual void defineBlockPC( std::vector<Preconditioner<GKVector,GKOps>*>&,
                                  std::vector<DOFList>&,
                                  const GKVector&,
                                  GKOps&,
                                  const std::string&,
                                  const std::string&,
                                  bool,
                                  const KineticSpecies&,
                                  const GlobalDOFKineticSpecies&,
                                  int ) { }

      virtual void updateBlockPC( std::vector<Preconditioner<GKVector,GKOps>*>&,
                                  const KineticSpecies&,
                                  const GlobalDOFKineticSpecies&,
                                  const Real,
                                  const Real,
                                  const bool,
                                  const int ) { }

      virtual void preTimeStep  (const KineticSpeciesPtrVect&,
                                 const int, const int, const Real,
                                 const KineticSpeciesPtrVect&) { }

      virtual void postTimeStage(const KineticSpeciesPtrVect& soln,
                                 const int, const int, const Real, const int ) { }
};

#include "NamespaceFooter.H"
#endif
