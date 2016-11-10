#ifndef CASM_ConfigEnumStrain
#define CASM_ConfigEnumStrain

#include "casm/symmetry/PermuteIterator.hh"
#include "casm/strain/StrainConverter.hh"
#include "casm/container/InputEnumerator.hh"
#include "casm/container/Counter.hh"
#include "casm/clex/Configuration.hh"

namespace CASM {


  ENUMERATOR_INTERFACE_TRAITS(ConfigEnumStrain)

  /// Enumerate strained Configurations
  ///
  /// \ingroup ConfigEnum
  ///
  class ConfigEnumStrain : public InputEnumeratorBase<Configuration> {

    // -- Required members -------------------

  public:

    ConfigEnumStrain(Supercell &scel,
                     const Configuration &_init,
                     const std::vector<Index> &subspace_partitions,
                     const std::vector<double> &magnitudes,
                     std::string _mode);

    ENUMERATOR_MEMBERS(ConfigEnumStrain)

  private:

    /// Implements increment over all strain states
    void increment() override;


    // -- Unique -------------------

    Configuration m_current;

    // counts over strain grid
    EigenCounter<Eigen::VectorXd> m_counter;
    // counts over transformation matrices
    Index m_equiv_ind;
    StrainConverter m_strain_calc;
    //set of non-equivalent transformation matrices matrices that, along with m_counter define irreducible space
    std::vector<Eigen::MatrixXd> m_trans_mats;
    PermuteIterator m_perm_begin, m_perm_end;
    Eigen::MatrixXd m_shape_factor;

    const PermuteIterator &_perm_begin() {
      return m_perm_begin;
    }
    const PermuteIterator &_perm_end() {
      return m_perm_end;
    }

  };

}

#endif
