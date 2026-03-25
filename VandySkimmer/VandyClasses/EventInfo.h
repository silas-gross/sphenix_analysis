#include <iostream>
#include <limits>
#include <vector>

class EventInfo {
    public:
        EventInfo() = default;
        virtual ~EventInfo() = default;

        double get_z_vtx() {return z_vtx;};
        void set_z_vtx(double vtx){ z_vtx = vtx; };

        double get_z_vtx_truth() {return z_vtx_truth;};
        void set_z_vtx_truth(double vtx){ z_vtx_truth = vtx; };        

        double get_ZDC_rate() {return ZDC_rate;};
        void set_ZDC_rate(double rate){ ZDC_rate = rate; };

        double get_cross_section() {return cross_section;};
        void set_cross_section(double cs){ cross_section = cs; };

        bool is_dijet_event(int jetR_index) { return dijet_event[jetR_index]; };
        void set_dijet_event(int jetR_index, bool dijet){ dijet_event[jetR_index] = dijet; };

        float get_lead_pT(int jetR_index) {return lead_pT[jetR_index];};
        void set_lead_pT(int jetR_index, float pT){ lead_pT[jetR_index] = pT; };

        float get_sublead_pT(int jetR_index) {return sublead_pT[jetR_index];};
        void set_sublead_pT(int jetR_index, float pT){ sublead_pT[jetR_index] = pT; };

        bool is_dijetTruth_event(int jetR_index) { return dijetTruth_event[jetR_index]; };
        void set_dijetTruth_event(int jetR_index, bool dijet){ dijetTruth_event[jetR_index] = dijet; };

        float get_leadTruth_pT(int jetR_index) {return leadTruth_pT[jetR_index];};
        void set_leadTruth_pT(int jetR_index, float pT){ leadTruth_pT[jetR_index] = pT; };

        float get_subleadTruth_pT(int jetR_index) {return subleadTruth_pT[jetR_index];};
        void set_subleadTruth_pT(int jetR_index, float pT){ subleadTruth_pT[jetR_index] = pT; };

    private:
        double z_vtx{0.0};
        double z_vtx_truth{0.0};
        double ZDC_rate{0.0};
        double cross_section{0.0};
        bool dijet_event[4] = {false, false, false, false};
        float lead_pT[4] = {-999, -999, -999, -999};
        float sublead_pT[4] = {-999, -999, -999, -999};
        bool dijetTruth_event[4] = {false, false, false, false};
        float leadTruth_pT[4] = {-999, -999, -999, -999};
        float subleadTruth_pT[4] = {-999, -999, -999, -999};

    };