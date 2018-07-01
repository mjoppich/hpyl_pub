import * as React from 'react'; 

import ChipInput from 'material-ui-chip-input'
import AutoComplete from 'material-ui/AutoComplete'
import axios from 'axios';
import config from '../config';

interface IEntityElement {
    name: string;
    id: string;
 }

export interface OrganismChipACProps { onValueChange: any};
export interface OrganismChipACState { values: Array<string>, datasource: Array<string>, selected_elements: Array<IEntityElement>, elements: Array<IEntityElement>};


export default class OrganismChipAC extends React.Component<OrganismChipACProps, OrganismChipACState>{

    neo4jd3: any = null;

    constructor(props)
    {
        super(props);
    }

    componentWillMount()
    {
        var self=this;
        this.setState({values: [], datasource: [], selected_elements:[]});

        axios.get(config.getRestAddress() + "/defaultorgs", config.axiosConfig)
        .then(function (response) {
            console.log(response.data)

            var allElements = response.data;
            var elemnames = [];
            allElements.forEach(element => {
                elemnames.push(element.name)
            });

            self.setState({values: elemnames, datasource: elemnames, elements: allElements, selected_elements: allElements})
            self.props.onValueChange(allElements);

        })
        .catch(function (error) {
            console.log(error)
            self.setState({elements: [], datasource: [], values: [], selected_elements: []})
            self.props.onValueChange([]);
        });
    }

    componentDidMount(){

    }

    acceptInput(value)
    {

        for (var i = 0; i < this.state.elements.length; ++i) {

            if (this.state.elements[i].name == value)
            {
                console.log("accepting input")
                console.log(value)
                return true;
            }
        }

        console.log("Not accepting input")
        console.log(value)
        return false;

    }

    addElement(value)
    {

        console.log("Going to add")
        console.log(value)

        for (var i = 0; i < this.state.elements.length; ++i) {

            if (this.state.elements[i].name == value)
            {
                var thisElem = this.state.elements[i];
                this.state.selected_elements.push(thisElem);
                this.state.values.push(this.state.elements[i].name);
                this.setState({values: this.state.values, selected_elements: this.state.selected_elements});
            }
        }

        console.log("ChipAC sending data");
        console.log(this.state.selected_elements);
        this.props.onValueChange(this.state.selected_elements);

    }
    deleteElement( elementText, i )
    {
        var idx = this.state.values.indexOf(elementText);

        if (idx >= 0)
        {
            this.state.values.splice(idx, 1);
            this.state.selected_elements.splice(idx, 1);
        }

        this.props.onValueChange(this.state.selected_elements);

    }

    handleOrganismAC(searchText)
    {
        var self = this;

        axios.post(config.getRestAddress() + "/orgautocomplete", {search: searchText}, config.axiosConfig)
        .then(function (response) {
            console.log(response.data)

            var allElements = response.data;
            var elemnames = [];
            allElements.forEach(element => {
                elemnames.push(element.name)
            });

            self.setState({datasource: elemnames, elements: allElements})

        })
        .catch(function (error) {
            console.log(error)
            self.setState({elements: [], datasource: []})
        });
    }
  
    render(){
      return (

                <ChipInput
                            value={this.state.values}
                            onBeforeRequestAdd={(chip) => this.acceptInput(chip)}
                            onRequestAdd={(chip) => this.addElement(chip)}
                            onRequestDelete={(chip, index) => this.deleteElement(chip, index)}
                            openOnFocus={true}
                            filter={AutoComplete.fuzzyFilter}
                            dataSource={this.state.datasource}
                            fullWidth
                            fullWidthInput
                            floatingLabelText='Selected Organisms'
                            hintText='Organism ID or Name'
                            onUpdateInput={(searchText, dataSource, params) => this.handleOrganismAC(searchText)}
                        />

        )
    }
}