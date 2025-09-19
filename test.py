import streamlit as st
import electromagnetism as em

st.title("My Python Package Interface")

# Add widgets to interact with your package
user_input = st.text_input("Enter input for your package:")
if st.button("Run"):
    result = em.models.coil.Coil(user_input)
    st.write("Result:", result)