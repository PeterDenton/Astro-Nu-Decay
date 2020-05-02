#ifndef gsl2dint_H
#define gsl2dint_H

// adapted from https://stackoverflow.com/questions/40157295/efficient-double-integration-in-c-with-gsl
// Simple RAII wrapper 
class IntegrationWorkspace
{
	gsl_integration_workspace *w;

	public:
	IntegrationWorkspace(const size_t n = 1000): w(gsl_integration_workspace_alloc(n)) {} 
	~IntegrationWorkspace()
	{
		gsl_integration_workspace_free(w);
	}

	operator gsl_integration_workspace*()
	{
		return w; 
	}
};

// Build gsl_function from lambda
template <typename F>
class gsl_function_pp: public gsl_function {
    const F func;
    static double invoke(double x, void *params) {
        return static_cast<gsl_function_pp*>(params)->func(x);
    }
    public:
    gsl_function_pp(const F& f) : func(f)
    {
        function = &gsl_function_pp::invoke; //inherited from gsl_function
        params   = this;                                         //inherited from gsl_function
    }
    operator gsl_function*(){return this;}
};

// Helper function for template construction
template <typename F> gsl_function_pp<F> make_gsl_function(const F& func)
{
    return gsl_function_pp<F>(func);
}

#endif
