# blinggen module 

This is a collection of modules designed to work with McStas inputs, parse them, 
and built an MCNP beamline model for shielding. It takes care of some rather 
tricky things, such as plane definition, mirror thickness (dependant on m-value),
and generates the transformates for different components of the guide, which is 
very useful for later edition, such as shielding creation and mess tallying. 

For now, only the McStas parser is working. Be aware that McStas is a pretty
wide metalanguage, so, when applied to other instruments, more functions will
probably need to be implemented. The MCNP input generator is a Work-In-Progress
now. For other codes, the same results from the parser can be used to generate
inputs

If you use this module and find it
useful, or have suggestions, issues, pull requests, etc, feel free to write.

CONTRIBUTING

Contributions are highly welcome! If you would like to contribute, it is as 
easy as forking the repository on GitHub, making your changes, and issuing 
a pull request. If you have any questions about this process don't hesitate 
to ask the author.
