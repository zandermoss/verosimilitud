class Flux
{
	
	public:

		Flux(double arg_normalization, double arg_index);
		double evaluatePowerlawFlux(const powerlawFlux flux, double neutrinoEnergy);
		double integratePowerlawFlux(const powerlawFlux flux, double minEnergy, double maxEnergy){
		void SetNormalization(double arg_normalization);
		void SetIndex(double arg_index);
		double GetNormalization(void);
		double GetIndex(void);

	private:

		double normalization;
		double index;
};
