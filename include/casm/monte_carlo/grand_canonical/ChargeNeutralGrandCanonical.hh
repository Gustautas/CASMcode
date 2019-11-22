#ifndef CASM_ChargeNeutralGrandCanonical_HH
#define CASM_ChargeNeutralGrandCanonical_HH

#include "casm/monte_carlo/MonteCarlo.hh"
#include "casm/monte_carlo/grand_canonical/GrandCanonicalConditions.hh"
#include "casm/monte_carlo/grand_canonical/GrandCanonicalSettings.hh"

namespace CASM {
  ///
  /// Derives from base MonteCarlo class, to be used for simulations at constant
  /// temperature and chemical potential and accounts for charge neutral swaps.
  ///
  /// As with all the other derived Monte Carlo classes, member functions must
  /// follow a specific naming convention to be used with templated routines currently
  /// defined in MonteDriver.hh:
  ///      -conditions (should be able to steal from GrandCanonical)
  ///      -set_conditions (should be able to steal from GrandCanonical)
  ///      -propose (we have to write this for our specialized event)
  ///      -check  (we have to write this for our specialized event)
  ///      -accept(we have to write this for our specialized event)
  ///      -reject(we have to write this for our specialized event)
  ///      -write_results
  ///
  //
  class ChargeNeutralGrandCanonicalEvent;


class ChargeNeutralGrandCanonical : public MonteCarlo{
	public:

    typedef ChargeNeutralGrandCanonicalEvent EventType;
    typedef GrandCanonicalConditions CondType;
    typedef GrandCanonicalSettings SettingsType;

	ChargeNeutralGrandCanonical(PrimClex &primclex, const SettingsType &settings, Log &_log);


    /// \brief Return current conditions
    const CondType &conditions() const;
    
    ///Keeps track of what sites can change to what
    const SiteExchanger m_site_swaps;

	/// \brief Set conditions and clear previously collected data
    void set_conditions(const CondType &new_conditions);

	/// This function needs to do all the math for energy and correlation deltas and store
	/// the results inside the containers hosted by event.
	void _update_deltas(EventType &event, 
						std::pair<Index,Index> &mutating_sites,
						std::pair<int,int> &sublats,
						std::pair<int,int> &curr_occs,
						std::pair<int,int> &new_occs) const;

    /// \brief Propose a new event, calculate delta properties, and return reference to it
    const EventType &propose();
    
	/// \brief Based on a random number, decide if the change in energy from the proposed event is low enough to be accepted.
    bool check(const EventType &event);
    
	/// \brief Accept proposed event. Change configuration accordingly and update energies etc.
    void accept(const EventType &event);

    /// \brief Nothing needs to be done to reject a GrandCanonicalEvent
    void reject(const EventType &event);
    
	/// \brief Write results to files
    void write_results(Index cond_index) const;

};




}
#endif